from iedb_receptor_pipeline.utilities import util
from iedb_receptor_pipeline.utilities.logging_setup import set_logging_context
import numpy as np
import logging
import pandas as pd


def safe_get_field(tool_data, field):
    result = tool_data[field] if tool_data is not None and field in tool_data else None

    if result in ("", "nan"):
        result = None
    elif isinstance(result, float) and np.isnan(result):
        result = None

    if type(result) == str:
        result = result.strip()

    return result

def safe_correct_gene(field):
    if field is None:
        return None

    return ", ".join(sorted([gene.strip() for gene in field.split(",")]))


def is_field_in_igblast_list(field_str, igblast_list_str):
    if igblast_list_str is not None and field_str is not None:
        igblast_genes = [util.remove_allele(igblast_field_str.strip()) for igblast_field_str in igblast_list_str.split(",")]
        return util.remove_allele(field_str) in igblast_genes
    return True


def get_nt_calc_data(input_row, igblast_nt_row):
    return {"receptor_type": safe_get_field(input_row, "receptor_type"),
            "chain_type": consolidate_chain_type(input_row, igblast_nt_row),
            "pro_seq": safe_get_field(input_row, "pro_seq_orig"),
            "v_gene": safe_correct_gene(safe_get_field(igblast_nt_row, "v_call")),
            "d_gene": safe_correct_gene(safe_get_field(igblast_nt_row, "d_call")),
            "j_gene": safe_correct_gene(safe_get_field(igblast_nt_row, "j_call")),
            "junction_seq": safe_get_field(igblast_nt_row, "junction_aa"),
            "cdr3_seq": safe_get_field(igblast_nt_row, "cdr3_aa"),
            "cdr2_seq": safe_get_field(igblast_nt_row, "cdr2_aa"),
            "cdr1_seq": safe_get_field(igblast_nt_row, "cdr1_aa"),
            "v_dom_seq": safe_get_field(igblast_nt_row, "sequence_alignment_aa"),
            "v_gene_method": "IgBLAST nt",
            "j_gene_method": "IgBLAST nt",
            "d_gene_method": "IgBLAST nt",
            "junction_seq_method":"IgBLAST nt",
            "cdr3_seq_method": "IgBLAST nt",
            "cdr2_seq_method": "IgBLAST nt",
            "cdr1_seq_method": "IgBLAST nt",
            "v_dom_seq_method": "IgBLAST nt"}


def consolidate_v_gene_aa(igblast_aa_row, tt_row):
    igblast_aa_v = safe_correct_gene(safe_get_field(igblast_aa_row, "v_call"))
    tt_v = safe_correct_gene(safe_get_field(tt_row, "vgene_calculated"))

    # if aa info is available, prefer IgBLAST aa over curated, warn only if gene-level is different
    if igblast_aa_v is not None:
        if tt_v is not None and igblast_aa_v != tt_v:
            if not is_field_in_igblast_list(tt_v, igblast_aa_v):
                logging.warning(
                    f"Mismatch between V gene calculated by IgBLAST amino acid ({igblast_aa_v}) "
                    f"and curated V gene corrected with tidytcells ({tt_v}). "
                    f"Calculated V gene will be determined by IgBLAST.")

        return igblast_aa_v, "IgBLAST aa"

    # only curated is available
    if tt_v is not None:
        return tt_v, "Curated, standardized"

    # no info available
    return None, None


def consolidate_cdr3_junction(input_row, anarcii_row, tt_row):
    input_seq = safe_get_field(input_row, "cdr3_seq_curated")
    anarcii_junction = safe_get_field(anarcii_row, "junction_seq_calculated")
    anarcii_cdr3 = safe_get_field(anarcii_row, "cdr3_seq_calculated")
    tt_junction = safe_get_field(tt_row, "junction_seq_calculated")
    tt_cdr3 = safe_get_field(tt_row, "cdr3_seq_calculated")

    if anarcii_cdr3:
        if tt_cdr3 is None and tt_junction is None:
            return anarcii_cdr3, anarcii_junction, "ANARCII", "ANARCII"
        if (tt_cdr3 is not None and anarcii_cdr3 != tt_cdr3) or (tt_junction is not None and anarcii_junction != tt_junction):
            logging.warning(f"Mismatch between CDR3/Junction by ANARCI ({anarcii_cdr3}, {anarcii_junction}) "
                          f"and curated CDR3/Junction corrected with tidytcells ({tt_cdr3}, {tt_junction}). "
                          F"The ANARCI CDR3/Junction will be used.")

            return anarcii_cdr3, anarcii_junction, "ANARCII", "ANARCII"
        elif (tt_cdr3 is not None and anarcii_cdr3 == tt_cdr3) and (tt_junction is not None and anarcii_junction == tt_junction):
            return anarcii_cdr3, anarcii_junction, "ANARCII", "ANARCII"

    if tt_cdr3:
        return tt_cdr3, tt_junction, "Curated, standardized", "Curated, standardized"

    if input_seq:
        if util.matches_basic_junction_pattern(input_seq):
            input_junction = util.safe_return_junction(input_seq, input_row["aa"])
            if input_junction:
                input_cdr3 = input_junction[1:-1]
                return input_cdr3, input_junction, "Curated", "Curated"

    return None, None, None, None


def consolidate_cdr12(input_row, anarcii_aa_row, tt_row, cdr):
    input_cdr = safe_get_field(input_row, f"{cdr}_seq_curated")
    anarcii_cdr = safe_get_field(anarcii_aa_row, f"{cdr}_seq_calculated")
    tt_cdr = safe_get_field(tt_row, f"{cdr}_seq_calculated")

    if anarcii_cdr is not None:
        if input_cdr is not None and anarcii_cdr != input_cdr:
            logging.warning(f"Mismatch between ANARCII {cdr.upper()} sequence {anarcii_cdr} and "
                          f"and curated {cdr.upper()} ({input_cdr}). "
                          f"The ANARCII {cdr.upper()} sequence will be used.")

        return anarcii_cdr, "ANARCII"

    if tt_cdr is not None:
        if input_cdr is None:
            return tt_cdr, "Germline"
        return tt_cdr, "Curated, standardized"

    if input_cdr is not None:
        input_cdr = input_cdr.strip().upper()
        input_cdr = input_cdr if all(aa in util.AA_ALPHABET_STRICT for aa in input_cdr) else None

        if input_cdr is not None:
            return input_cdr, "Curated"

    return None, None


def check_gene_type(gene, gene_method, chain_type, gene_type):
    if gene is not None and chain_type is not None:
        if len(gene) < 4 or gene[0:3] not in util.CHAIN_TO_LEGAL_LOCI[chain_type] or gene[3] != gene_type:
            if chain_type == "delta" and gene[0:4] == "TRAV" and "DV" in gene:
                pass
            else:
                logging.error(f"Consolidated {chain_type} chain {gene_type} ({gene_method}) gene seems incorrect: '{gene}' (expected '{'/'.join(list(util.CHAIN_TO_LEGAL_LOCI[chain_type]))}{gene_type}...')")


def consolidate_chain_type(input_row, igblast_row, anarcii_aa_row=None, v_gene=None):
    if igblast_row is not None:
        chain_type = util.CHAIN_LOCUS_TO_IEDB[igblast_row["locus"]]
        method = "IgBLAST"
    elif anarcii_aa_row is not None:
        chain_type = safe_get_field(anarcii_aa_row, "chain_type")
        method = "ANARCII"
    else:
        chain_type = safe_get_field(input_row, "chain_type")
        method = None

    if chain_type != input_row["chain_type"]:
        if not chain_type in ('kappa_light', 'lambda_light') and input_row["chain_type"] == 'light':
            logging.error(f"{method} determined the chain type to be {chain_type} instead of curator-assigned chain type {input_row['chain_type']}. The chain type has been corrected.")

    if chain_type == "light" and v_gene is not None:
        chain_type = util.CHAIN_LOCUS_TO_IEDB[v_gene[0:3]]

    valid_chains_for_receptor = util.RECEPTOR_TYPE_TO_VALID_CHAINS[str(input_row["receptor_type"]).lower()]

    if chain_type is not None and chain_type.lower() not in valid_chains_for_receptor:
        logging.error(f"Chain type {chain_type} is not recognized for receptor type {input_row['receptor_type']}.")

    return chain_type

def get_aa_calc_data(input_row, igblast_aa_row, anarcii_aa_row, tt_row):
    v_gene, v_gene_method = consolidate_v_gene_aa(igblast_aa_row, tt_row)
    d_gene = safe_correct_gene(safe_get_field(tt_row, "dgene_calculated"))
    d_gene_method = "Curated, standardized" if d_gene is not None else None
    j_gene = safe_correct_gene(safe_get_field(tt_row, "jgene_calculated"))
    j_gene_method = "Curated, standardized" if j_gene is not None else None

    chain_type = consolidate_chain_type(input_row, igblast_aa_row, anarcii_aa_row=anarcii_aa_row, v_gene=v_gene)

    check_gene_type(v_gene, v_gene_method, chain_type=chain_type, gene_type="V")
    check_gene_type(d_gene, d_gene_method, chain_type=chain_type, gene_type="D")
    check_gene_type(j_gene, j_gene_method, chain_type=chain_type, gene_type="J")

    cdr3, junction, cdr3_method, junction_method = consolidate_cdr3_junction(input_row, anarcii_aa_row, tt_row)
    cdr1, cdr1_method = consolidate_cdr12(input_row, anarcii_aa_row, tt_row, "cdr1")
    cdr2, cdr2_method = consolidate_cdr12(input_row, anarcii_aa_row, tt_row, "cdr2")
    v_dom_seq = safe_get_field(anarcii_aa_row, "v_dom_seq")
    v_dom_seq_method = "ANARCII" if v_dom_seq is not None else None

    return {"receptor_type": safe_get_field(input_row, "receptor_type"),
            "chain_type": chain_type,
            "pro_seq": safe_get_field(input_row, "pro_seq_orig"),
            "v_gene": v_gene, "v_gene_method": v_gene_method,
            "d_gene": d_gene, "d_gene_method": d_gene_method,
            "j_gene": j_gene, "j_gene_method": j_gene_method,
            "junction_seq": junction, "junction_seq_method": junction_method,
            "cdr3_seq": cdr3, "cdr3_seq_method": cdr3_method,
            "cdr2_seq": cdr2, "cdr2_seq_method": cdr2_method,
            "cdr1_seq": cdr1, "cdr1_seq_method": cdr1_method,
            "v_dom_seq": v_dom_seq, "v_dom_seq_method": v_dom_seq_method}


def to_indexed(df):
    if df is not None:
        return df.set_index(["row", "chain"], drop=False) # if {"row","chain"} <= set(df.columns) else df

def get_row(df, row, chain):
    if df is not None:
        try:
            return df.loc[(row, chain)]
        except KeyError:
            return None


def consolidate_results(tool_input, anarcii_results, igblast_nt_results, igblast_aa_results, tt_results):
    tool_input = to_indexed(tool_input)
    tool_input["species_ncbi"] = tool_input["species_latin"].map(util.LATIN_TO_NCBITAXON)

    anarcii_results = to_indexed(anarcii_results)
    igblast_nt_results = to_indexed(igblast_nt_results)
    igblast_aa_results = to_indexed(igblast_aa_results)
    tt_results = to_indexed(tt_results)

    df_rows = []
    total = len(tool_input)

    for i, row in enumerate(tool_input["row"].unique()):
        if i % 100 == 0:
            print(f"Consolidating {round(i/total*100,1)}%")

        for chain in (1, 2):
            with set_logging_context(template_row=row, chain=chain):
                input_row = get_row(tool_input, row, chain)
                if input_row is None:
                    continue

                species = input_row["species_ncbi"]

                igblast_nt_row = get_row(igblast_nt_results, row, chain)
                if igblast_nt_row is not None:
                    df_rows.append({"row": row, "chain": chain, "species": species,
                                    **get_nt_calc_data(input_row, igblast_nt_row)})
                    continue

                igblast_aa_row = get_row(igblast_aa_results, row, chain)
                anarcii_row = get_row(anarcii_results, row, chain)
                tt_row = get_row(tt_results, row, chain)

                df_rows.append({"row": row, "chain": chain, "species": species,
                                **get_aa_calc_data(input_row,
                                                   igblast_aa_row,
                                                   anarcii_row,
                                                   tt_row)})


    print(f"Consolidating 100.0%")

    results = pd.DataFrame(df_rows)

    return results
