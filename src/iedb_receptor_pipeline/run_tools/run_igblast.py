import subprocess

from iedb_receptor_pipeline.utilities.logging_setup import set_logging_context
import iedb_receptor_pipeline.utilities.util as util
import logging
from pathlib import Path
import pandas as pd


def row_to_fasta(row, seq_col="nt"):
    return f">{row['row']}_{row['chain']}\n{row[seq_col]}\n"


def chain_to_igblast_receptor(chain):
    if chain in util.TCR_CHAINS:
        return "TCR"
    elif chain in util.BCR_CHAINS:
        return "Ig"
    else:
        return "NA"

def nt_is_cdr3(row):
    run_igblast = True

    with set_logging_context(template_row=row["row"], chain=row["chain"]):
        if len(row["nt"]) <= 120:
            run_igblast = False

            if not pd.isna(row["cdr3_seq_curated"]) and len(row["cdr3_seq_curated"]) * 3 == len(row["nt"]):
                pass
            else:
                found = "NA" if pd.isna(row["cdr3_seq_curated"]) else row["cdr3_seq_curated"]
                logging.info(
                    f"Skipping IgBLAST nt because the input nucleotide sequence {row["nt"]} is too short (<= 120 nucleotides). If this is the CDR3 sequence, please add the translated sequence as the Curated CDR3 sequence instead (found {found}).")

    return run_igblast

def aa_is_cdr3(row):
    run_igblast = True

    with set_logging_context(template_row=row["row"], chain=row["chain"]):
        if len(row["aa"]) <= 40:
            run_igblast = False

            if row["aa"] == row["cdr3_seq_curated"]:
                logging.info(f"The full protein sequence for this chain is identical to the curated CDR3 sequence. It is preferred to only keep the curated CDR3 sequence.")
            else:
                found = "NA" if pd.isna(row["cdr3_seq_curated"]) else row["cdr3_seq_curated"]
                logging.info(f"Skipping IgBLAST aa because the input protein sequence {row["aa"]} is too short (<= 40 amino acids). If this is the CDR3 sequence, please add this as the Curated CDR3 sequence instead (found {found}).")

    return run_igblast

def get_igblast_input(full_input, seq_col):
    igblast_input = full_input[full_input[seq_col].notna()][["row", "chain", "species_latin", seq_col, "cdr3_seq_curated", "chain_type", "receptor_type"]].copy()
    igblast_input["species_igblast"] = igblast_input["species_latin"].map(lambda x: util.LATIN_TO_COMMON[x] if x in util.LATIN_TO_COMMON else "NA")
    igblast_input["receptor_igblast"] = igblast_input["chain_type"].map(chain_to_igblast_receptor)

    igblast_input = igblast_input.loc[(igblast_input["species_igblast"] != "NA") &
                                      (igblast_input["receptor_igblast"] != "NA") &
                                      (igblast_input["chain_type"] != "construct")]

    if seq_col == "aa":
        run_igblast = igblast_input.apply(aa_is_cdr3, axis=1)
    elif seq_col == "nt":
        run_igblast = igblast_input.apply(nt_is_cdr3, axis=1)

    return igblast_input[run_igblast].copy()

def run_query(query):
    logging.error("No query provided")
    if query is None:
        return False

    logging.debug("Running query:\n" + query)
    result = subprocess.run(query, shell=True, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logging.error(result.stderr)
        return False

    return True


def get_db_name(receptor_type, species):
    if receptor_type == "TCR" or species in ("mouse", "rhesus_monkey"):
        return "imgt"

    if receptor_type == "Ig" and species == "human":
        return "ogrdb"

    assert False, "Unknown db"


def construct_igblast_nt_query(species, receptor_type, igblast_input, igblast_output, tmp_folder, db_name, n_processes=2):
    assert species in {"human", "mouse", "rhesus_monkey"}
    assert receptor_type in {"Ig", "TCR"}
    rt = "tr" if receptor_type == "TCR" else "ig"

    GERMLINE_ROOT = util.OGRDB_GERMLINE_ROOT if db_name == "ogrdb" else util.IMGT_GERMLINE_ROOT

    if GERMLINE_ROOT is None:
        return None

    germline_folder = GERMLINE_ROOT / "database"
    optional_folder = GERMLINE_ROOT / "optional_file"

    # note: always use imgt C region because we dont have them for OGRDB
    query = f'''export IGDATA={util.IGBLAST_ROOT} && {util.IGBLAST_N} \
              -germline_db_V {germline_folder}/{db_name}_{species}_{rt}_v \
              -germline_db_D {germline_folder}/{db_name}_{species}_{rt}_d \
              -germline_db_J {germline_folder}/{db_name}_{species}_{rt}_j \
              -c_region_db {util.IMGT_GERMLINE_ROOT / 'database'}/imgt_{species}_{rt}_c \
              -auxiliary_data {optional_folder}/{species}_gl.aux \
              -domain_system imgt -ig_seqtype {receptor_type} -organism {species} \
              -outfmt 19 -num_threads {n_processes} \
              -query {tmp_folder}/{igblast_input.name} \
              -out {tmp_folder}/{igblast_output.name}'''

    return query


def safe_merge_airr_results_and_tool_input(airr_df, input_df):
    column_overlap = set(input_df.columns).intersection(set(airr_df.columns))

    assert len(column_overlap) == 0, \
        f"Cannot safely merge IgBLAST input and AIRR output due to overlapping columns: {column_overlap}"
    assert "row" not in airr_df and "chain" not in airr_df, \
        f"Cannot safely merge IgBLAST input and AIRR output due to AIRR file containing columns named 'row' or 'chain'."

    airr_df[["row", "chain"]] = airr_df["sequence_id"].str.split("_", expand=True)
    airr_df["row"] = airr_df["row"].astype(str)
    airr_df["chain"] = airr_df["chain"].astype(int)

    return input_df.merge(airr_df, on=["row", "chain"])

def re_compute_locus(row):
    igblast_locus = util.safe_get_field_from_df(row["locus"])

    if igblast_locus in ("TRA", "TRD") and row["species_latin"] not in ("Homo sapiens", "Mus musculus", "Macaca mulatta"):
        logging.warning(f"Recompute locus: Not sure how to deal with species {row["species_latin"]}, returning IgBLAST locus {igblast_locus}")
        return igblast_locus

    v_call = util.safe_get_field_from_df(row["v_call"])
    d_call = util.safe_get_field_from_df(row["d_call"])
    j_call = util.safe_get_field_from_df(row["j_call"])
    curator_assigned = util.safe_get_field_from_df(row["chain_type"])
    receptor_type = util.safe_get_field_from_df(row["receptor_type"])

    if igblast_locus == "TRA" and curator_assigned == "delta" and receptor_type == "gammadelta":
        if ((v_call is not None and v_call.startswith("TRA") and "/DV" in v_call) or
                (d_call is not None and d_call.startswith("TRD")) or
                (j_call is not None and j_call.startswith("TRD"))):
            logging.info("IgBLAST/MakeDb.py assigned TRA locus to a chain curated as 'delta' in an 'gammadelta' receptor due to shared TRAV/DV gene. Correcting IgBLAST locus to 'TRD'.")
            return "TRD"

    if igblast_locus == "TRD" and curator_assigned == "alpha" and receptor_type == "alphabeta":
        logging.info("IgBLAST/MakeDb.py assigned TRD locus to a chain curated as 'alpha' in an 'alphabeta' receptor. TRD genes are commonly known to be used by TRA chains. Correcting IgBLAST locus to 'TRA'.")
        return "TRA"

    return igblast_locus


def remove_poor_quality(igblast_output, iblast_type):
    filtered_rows = []
    for index, curr_row in igblast_output.iterrows():
        with set_logging_context(template_row=curr_row["row"], chain=curr_row["chain"]):
            if pd.isna(curr_row["locus"]):
                logging.error(f"Failed to compute IgBLAST {iblast_type} results for this chain. This may indicate the input sequence is not a receptor chain.")
                continue

            curr_row["locus"] = re_compute_locus(curr_row)

            if curr_row["locus"] not in util.CHAIN_TO_LEGAL_LOCI[curr_row["chain_type"]]:
                logging.error(f"Locus determined by IgBLAST {iblast_type} ({curr_row['locus']}, mapped to {util.CHAIN_LOCUS_TO_IEDB[curr_row['locus']]}) "
                                f"does not match the curator-assigned chain type ({curr_row['chain_type']}).")

            if "productive" in curr_row and curr_row["productive"] == "F":
                logging.info(f"IgBLAST {iblast_type} determined this rearrangement to be unproductive (stop_codon={curr_row['stop_codon']}, vj_in_frame={curr_row['vj_in_frame']}, v_frameshift={curr_row['v_frameshift']}).")

            if curr_row["v_support"] > util.EVAL_TRUST_CUTOFF:
                logging.warning(
                    f"The IgBLAST {iblast_type} alignment had V gene e-value {curr_row["v_support"]} above the trusted cutoff ({util.EVAL_TRUST_CUTOFF}). "
                    f"This may indicate that the curated species ({curr_row["species_igblast"]}) or chain type ({curr_row["chain_type"]}) were incorrectly assigned. "
                    f"The IgBLAST {iblast_type} V gene results for this chain will be skipped.")
                curr_row["vgene_calculated"] = None

            if "j_support" in curr_row and curr_row["j_support"] > util.EVAL_TRUST_CUTOFF:
                logging.warning(
                    f"The IgBLAST {iblast_type} alignment had J gene e-value {curr_row["j_support"]} above the trusted cutoff ({util.EVAL_TRUST_CUTOFF}). "
                    f"This may indicate that the curated species ({curr_row["species_igblast"]}) or chain type ({curr_row["chain_type"]}) were incorrectly assigned. "
                    f"The IgBLAST {iblast_type} J gene results for this chain will be skipped.")
                curr_row["jgene_calculated"] = None

            filtered_rows.append(curr_row)

    return pd.DataFrame(filtered_rows, columns=igblast_output.columns)


def process_igblast_results(airr_results, igblast_input, igblast_type):
    results = remove_poor_quality(safe_merge_airr_results_and_tool_input(airr_results, igblast_input), igblast_type)

    return results


def run_igblast_nt(full_input, tmp_folder, n_processes):
    igblast_nt_input = get_igblast_input(full_input, seq_col="nt")

    airr_results = []

    for (species, receptor_type), subset_input in igblast_nt_input.groupby(['species_igblast', 'receptor_igblast']):
        if species == "NA":
            logging.info(f"Skipping IgBLAST nucleotide for rows {util.str_format_unique_list(subset_input['row'])} due to unsupported species {util.str_format_unique_list(subset_input['species_latin'])}.")
            continue

        igblast_input = Path(tmp_folder) / f"igblast_nt_{species}_{receptor_type}.fasta"
        igblast_airr_output = Path(tmp_folder) / f"igblast_nt_{species}_{receptor_type}_airr.tsv"
        db_name = get_db_name(receptor_type, species)

        with open(igblast_input, "w") as f:
            f.writelines(list(subset_input.apply(row_to_fasta, axis=1, seq_col="nt")))

        query = construct_igblast_nt_query(species, receptor_type, igblast_input, igblast_airr_output, tmp_folder, db_name, n_processes)

        with set_logging_context(log_context="IgBLAST nucleotide"):
            igblast_success = run_query(query)

        if igblast_success:
            airr_results.append(pd.read_csv(igblast_airr_output, sep="\t"))

    if len(airr_results) >0:
        return process_igblast_results(pd.concat(airr_results), igblast_nt_input, "nucleotide")

def construct_igblast_aa_query(species, receptor_type, igblast_input, igblast_output, tmp_folder, db_name, n_processes=2):
    assert species in {"human", "mouse", "rhesus_monkey"}
    assert receptor_type in {"Ig", "TCR"}
    rt = "tr" if receptor_type == "TCR" else "ig"

    GERMLINE_ROOT = util.OGRDB_GERMLINE_ROOT if db_name == "ogrdb" else util.IMGT_GERMLINE_ROOT

    if GERMLINE_ROOT is None:
        return None

    germline_folder = GERMLINE_ROOT / "database"

    query = f'''export IGDATA={util.IGBLAST_ROOT} && {util.IGBLAST_P} \
                 -germline_db_V {germline_folder}/{db_name}_aa_{species}_{rt}_v \
                 -domain_system imgt -ig_seqtype {receptor_type} -organism {species} \
                 -outfmt '7 std qseq sseq btop' \
                 -num_threads {n_processes} \
                 -query {tmp_folder}/{igblast_input.name} \
                 -out {tmp_folder}/{igblast_output.name}'''

    return query

def construct_makedb_query(igblast_input, igblast_output, tmp_folder, airr_prefix, db_name):
    GERMLINE_ROOT = util.OGRDB_GERMLINE_ROOT if db_name == "ogrdb" else util.IMGT_GERMLINE_ROOT

    return f'''MakeDb.py igblast-aa --outdir {tmp_folder} \
                    --failed --extended \
                    --outname {airr_prefix} \
                    -s {tmp_folder}/{igblast_input.name} \
                    -i {tmp_folder}/{igblast_output.name} \
                    -r {GERMLINE_ROOT/ "vdj_aa"} \
                    --partial --extended'''


def run_igblast_aa(tool_input, tmp_folder, n_processes):
    igblast_aa_input = get_igblast_input(tool_input, seq_col="aa")

    airr_results = []

    for (species, receptor_type), subset_input in igblast_aa_input.groupby(['species_igblast', 'receptor_igblast']):
        if species == "NA":
            logging.info(f"Skipping IgBLAST protein for rows {util.str_format_unique_list(subset_input['row'])} due to unsupported species {util.str_format_unique_list(subset_input['species_latin'])}.")
            continue

        base_name = f"igblast_aa_{species}_{receptor_type}"
        igblast_input = Path(tmp_folder) / f"{base_name}.fasta"
        igblast_output = Path(tmp_folder) / f"{base_name}.txt"

        db_name = get_db_name(receptor_type, species)

        with open(igblast_input, "w") as f:
            f.writelines(list(subset_input.apply(row_to_fasta, axis=1, seq_col="aa")))

        query = construct_igblast_aa_query(species, receptor_type, igblast_input, igblast_output, tmp_folder, db_name, n_processes)

        with set_logging_context(log_context="IgBLAST protein"):
            igblast_success = run_query(query)

        if igblast_success:
            query = construct_makedb_query(igblast_input, igblast_output, tmp_folder, base_name + "_airr", db_name)

            with set_logging_context(log_context="MakeDb.py"):
                makedb_success = run_query(query)

            if makedb_success:
                for airr_file in Path(tmp_folder).glob(base_name + "_airr*"):
                    airr_results.append(pd.read_csv(airr_file,  sep="\t"))

    if len(airr_results) > 0:
        return process_igblast_results(pd.concat(airr_results), igblast_aa_input, "protein")

def test_igblast_setup():
    if util.IGBLAST_ROOT is None:
        mssg = "IgBLAST is not set up, please make sure util.IGBLAST_ROOT points to a valid IgBLAST source folder (e.g., 'ncbi-igblast-1.22.0')"
        logging.error(mssg)
        print(mssg)
        return False

    germline_mssg = "{} germline root folder is not set up, please make sure util.{}_GERMLINE_ROOT points to a folder containing the germline files with the following structure:\n" \
                     "germline_root_folder/\n" \
                     "  database/\n" \
                     "    ... IgBLAST database files \n" \
                     "  optional_file/\n" \
                     "    ... IgBLAST optional files"

    if util.OGRDB_GERMLINE_ROOT is None:
        mssg = germline_mssg.format("OGRDB", "OGRDB")
        logging.error(mssg)
        print(mssg)
        return False

    if util.IMGT_GERMLINE_ROOT is None:
        mssg = germline_mssg.format("IMGT", "IMGT") + "\n\nThe IMGT germline root is optional as long as the OGRDB root is set. Note that OGRDB is only set up for human IG germline reference sets at this time."
        logging.error(mssg)
        print(mssg)

    return True


def run_igblast(tool_input, tmp_folder, n_processes):
    if test_igblast_setup():
        nt_output = run_igblast_nt(tool_input, tmp_folder, n_processes)
        aa_output = run_igblast_aa(tool_input, tmp_folder, n_processes)

        return nt_output, aa_output

    return None, None
