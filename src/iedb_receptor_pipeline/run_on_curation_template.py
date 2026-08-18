from argparse import ArgumentParser

from parse_input import get_tool_input
from iedb_receptor_pipeline.run_tools.run_tools import run, setup
from iedb_receptor_pipeline.utilities.logging_setup import log_with_error_collector

import numpy as np
import pandas as pd
from pathlib import Path
from iedb_receptor_pipeline.utilities.logging_setup import set_logging_context, LOG_FILE_SEPARATOR
from iedb_receptor_pipeline.utilities import util
import logging

from iedb_receptor_pipeline.utilities.util import CURATION_TEMPLATE_COLUMNS

LOG_TEMPLATE_ROW, LOG_TEMPLATE_CHAIN = None, None

def parse_arguments():
    args = ArgumentParser()

    args.add_argument(f'--input_file', required=True)
    args.add_argument(f'--output_folder',  default="./iedb_receptors_output")

    return args.parse_args()


def safe_read_excel(input_file):
    full_input_orig = pd.read_excel(input_file, engine='openpyxl', dtype=str)
    full_input_orig = full_input_orig.replace('\ufeff', '', regex=True)

    return full_input_orig


def get_accession_to_seq_dict(accession_ids, seq_type):
    assert seq_type in ("nucleotide", "protein")

    if len(accession_ids) == 0:
        return dict()

    from Bio import Entrez # only import Entrez if needed, limit requests
    Entrez.email = None    # personal email may be set

    seq_dict = dict()

    with Entrez.efetch(id=accession_ids, db=seq_type, rettype="gb", retmode="xml") as handle:
        for record in Entrez.parse(handle):
            if 'GBSeq_sequence' in record:
                accession = record['GBSeq_accession-version'] if record['GBSeq_accession-version'] in accession_ids else record['GBSeq_primary-accession']
                seq_dict[accession] = record['GBSeq_sequence']

    return seq_dict


def get_accession_to_species(accession_ids, field='GBSeq_organism'):
    if len(accession_ids) == 0:
        return dict()

    from Bio import Entrez # only import Entrez if needed, limit requests
    Entrez.email = None    # personal email may be set

    with Entrez.efetch(id=accession_ids, db="protein", rettype="gb", retmode="xml") as handle:
        return {id: record[field] for id, record in zip(accession_ids, Entrez.parse(handle)) if 'GBSeq_sequence' in record}

def resolve_nt_for_row(curator_nt_seq):
    if not pd.isna(curator_nt_seq):
        return curator_nt_seq.lower().replace(" ", "")

def safe_get_aa_sequence(accession, accession_seq_map):
    if pd.isna(accession):
        return None

    if accession not in accession_seq_map:
        logging.error(f"No amino acid accession sequence was retrieved for {accession}")
        return None

    seq = str(accession_seq_map[accession]).upper().strip()

    if set(seq) == {"X"} or not util.is_valid_alphabet(seq, util.AA_ALPHABET_AMBIGUOUS):
        logging.error(f"Amino acid accession sequence for {accession} is not valid: {seq}")
        return None

    return seq

def safe_get_nt_sequence(accession, accession_seq_map):
    if pd.isna(accession):
        return None

    if accession not in accession_seq_map:
        logging.error(f"No nucleotide accession sequence was retrieved for {accession}")
        return None

    seq = str(accession_seq_map[accession]).lower().strip()

    if set(seq) == {"n"} or not util.is_valid_alphabet(seq, util.NT_ALPHABET_AMBIGUOUS):
        logging.error(f"Nucleotide accession sequence for {accession} is not valid: {seq}")
        return None

    return seq


def get_accession_chain_seq(row, chain, accession_seq_map, seq_type):
    if seq_type == "protein":
        seq_col = f"chain{chain}_pro_seq"
        acc_col = f"chain{chain}_accession"
        safe_get_sequence = safe_get_aa_sequence
    else:
        seq_col = f"chain{chain}_nucleotide"
        acc_col = f"chain{chain}_nt_accession"
        safe_get_sequence = safe_get_nt_sequence

    with set_logging_context(template_row=row.name, chain=chain):
        if pd.notna(row[seq_col]):
            curator_seq = row[seq_col].upper() if seq_type == "protein" else row[seq_col].lower()
        else:
            curator_seq = None

        # curator_seq = row[seq_col].upper() if pd.notna(row[seq_col]) else None

        accession = row[acc_col]
        accession_seq = safe_get_sequence(accession, accession_seq_map)

        if pd.notna(accession_seq):
            if pd.notna(curator_seq) and curator_seq != accession_seq:
                logging.warning( f"{seq_type.title()} sequence with accession {accession} does not match curator input (continuing with accession sequence)."
                                 f" Database sequence: {accession_seq}{LOG_FILE_SEPARATOR} Curator sequence: {curator_seq}")

            if curator_seq is None:
                logging.info(f"Curated {seq_type} sequence was retrieved from accession {accession}.")

            return accession_seq

        return curator_seq


def resolve_accession_sequences(curation_template, retrieve_accessions=True):
    if not retrieve_accessions:
        return None

    accession_ids_nt = {*curation_template[f"chain1_nt_accession"].dropna().unique(),
                     *curation_template[f"chain2_nt_accession"].dropna().unique()}

    accession_ids_pro = {*curation_template[f"chain1_accession"].dropna().unique(),
                     *curation_template[f"chain2_accession"].dropna().unique()}

    accession_seq_map_nt = get_accession_to_seq_dict(accession_ids_nt, "nucleotide")
    accession_seq_map_pro = get_accession_to_seq_dict(accession_ids_pro, "protein")

    for chain in (1, 2):
        curation_template[f"chain{chain}_nucleotide"] = curation_template.apply(lambda row: get_accession_chain_seq(row, chain, accession_seq_map_nt, "nucleotide"), axis=1)
        curation_template[f"chain{chain}_pro_seq"] = curation_template.apply(lambda row: get_accession_chain_seq(row, chain, accession_seq_map_pro, "protein"), axis=1)

def get_calculated_receptors_from_results(consolidate_results_df):
    consolidate_results_df["row"] = consolidate_results_df["row"].astype("str")
    calculated_receptors = consolidate_results_df.pivot(index="row", columns="chain")
    calculated_receptors.columns = [f"{col}_chain{chain}" for col, chain in calculated_receptors.columns]
    calculated_receptors.reset_index(inplace=True)

    calculated_receptors["row"] = calculated_receptors["row"].astype(int)
    calculated_receptors = calculated_receptors.set_index("row").sort_index()
    calculated_receptors = calculated_receptors.replace({None: np.nan}).infer_objects(copy=False)

    return calculated_receptors

def remove_calculated_data_from_curation_template(curation_template):
    for field in util.CURATION_TEMPLATE_CALCULATED_COLUMNS:
        curation_template[field] = np.nan

def remove_obsolete_fields_curation_template(curation_template, drop_fields=True):
    '''This function safely drops obsolete columns known to occur in the old version of the curation template'''

    if not drop_fields:
        return

    cols_to_drop = [field for field in list(curation_template.columns)
                    if ("_start_" in field or "_end_" in field or field in ("Unnamed: 19", "calculated_receptor_type"))
                    and field not in CURATION_TEMPLATE_COLUMNS]

    for field in cols_to_drop:
        if not all(curation_template[field].isna()):
            logging.error(f"attempting to drop {field} but this is not NA!")
            print(f"ERROR: attempting to drop {field} but this is not NA!")
            drop_fields = False

    if drop_fields:
        logging.info(f"dropping the following fields: {cols_to_drop}")
        curation_template.drop(columns=cols_to_drop, inplace=True)

def rename_vdomain_field_to_curated(curation_template, rename=True):
    if not rename:
        return

    for chain in ("1", "2"):
        if f"chain{chain}_vdomain" in curation_template.columns:
            logging.info(f"Renaming colun 'chain{chain}_vdomain' to 'chain{chain}_vdomain_curated'")
            curation_template.rename(columns={f"chain{chain}_vdomain": f"chain{chain}_vdomain_curated"}, inplace=True)

def add_mandatory_fields_to_curation_template(curation_template):
    for col in util.CURATION_TEMPLATE_CURATED_COLUMNS:
        if col not in curation_template.columns:
            print(f"Adding missing column {col}")
            curation_template[col] = np.nan

def add_calculated_data_to_curation_template(curation_template, calculated_receptors):
    for chain in ("1", "2"):
        curation_template[f"chain{chain}_pro_seq"] = calculated_receptors[f"pro_seq_chain{chain}"]
        curation_template[f"chain{chain}_v_gene_calculated"] = calculated_receptors[f"v_gene_chain{chain}"]
        curation_template[f"chain{chain}_d_gene_calculated"] = calculated_receptors[f"d_gene_chain{chain}"]
        curation_template[f"chain{chain}_j_gene_calculated"] = calculated_receptors[f"j_gene_chain{chain}"]
        curation_template[f"chain{chain}_cdr1_seq_calculated"] = calculated_receptors[f"cdr1_seq_chain{chain}"]
        curation_template[f"chain{chain}_cdr2_seq_calculated"] = calculated_receptors[f"cdr2_seq_chain{chain}"]
        curation_template[f"chain{chain}_cdr3_seq_calculated"] = calculated_receptors[f"cdr3_seq_chain{chain}"]
        curation_template[f"chain{chain}_junction_seq_calculated"] = calculated_receptors[f"junction_seq_chain{chain}"]
        curation_template[f"chain{chain}_vdomain_calculated"] = calculated_receptors[f"v_dom_seq_chain{chain}"]
        curation_template[f"chain{chain}_chain_type_calculated"] = calculated_receptors[f"chain_type_chain{chain}"]

def add_errors_to_curation_template(curation_template, error_collector):
    curation_template["CHAIN1_ERRORS"] = curation_template.index.map(error_collector.errors_per_row["CHAIN1"])
    curation_template["CHAIN2_ERRORS"] = curation_template.index.map(error_collector.errors_per_row["CHAIN2"])


def reorder_columns(curation_template):
    extra_cols = [c for c in curation_template.columns if c not in util.CURATION_TEMPLATE_COLUMNS]
    curation_template = curation_template[util.CURATION_TEMPLATE_COLUMNS + extra_cols]

    return curation_template


def write_output(curation_template, output_folder):
    curation_template.to_excel(output_folder / f"{Path(args.input_file).stem}_calc.xlsx", index=False)

    erroneous_only = curation_template[curation_template["CHAIN1_ERRORS"].notna() | curation_template["CHAIN2_ERRORS"].notna()]
    erroneous_only.to_excel(output_folder / f"{Path(args.input_file).stem}_err_only.xlsx", index=False)

def main(args):
    output_folder, tmp_folder, intermediate_folder = setup(args.output_folder)

    curation_template = safe_read_excel(args.input_file)

    tmp_receptor_id = False
    if all(curation_template.curated_receptor_id.isna()):
        tmp_receptor_id = True
        curation_template.curated_receptor_id = curation_template.index.astype(str)

    remove_calculated_data_from_curation_template(curation_template)
    remove_obsolete_fields_curation_template(curation_template)
    rename_vdomain_field_to_curated(curation_template)
    add_mandatory_fields_to_curation_template(curation_template)

    with log_with_error_collector() as error_collector:
        resolve_accession_sequences(curation_template, retrieve_accessions=True)

        tool_input = get_tool_input(curation_template, translate_short_nt_to_cdr3_aa=True, v_dom_to_aa=True, scfv_dom_to_aa=True)


        consolidate_results_df = run(intermediate_folder, tmp_folder, tool_input, compute_results=True)
        calculated_receptors = get_calculated_receptors_from_results(consolidate_results_df)
        add_calculated_data_to_curation_template(curation_template, calculated_receptors)
        add_errors_to_curation_template(curation_template, error_collector)

    if tmp_receptor_id:
        curation_template.curated_receptor_id = np.nan

    curation_template = reorder_columns(curation_template)
    write_output(curation_template, output_folder)


if __name__ == '__main__':
    args = parse_arguments()
    main(args)
