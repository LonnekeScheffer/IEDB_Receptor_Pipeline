from argparse import ArgumentParser
from iedb_receptor_pipeline.run_tools.run_tools import setup, run
from iedb_receptor_pipeline.utilities import util
from iedb_receptor_pipeline.utilities.logging_setup import set_logging_context

from parse_input import format_tool_input_columns, validate_row_db

import pandas as pd
import logging

def parse_arguments():
    args = ArgumentParser()
    args.add_argument(f'--output_folder', default="./iedb_receptors_output")
    args.add_argument(f'--input_file_chain1', required=True)
    args.add_argument(f'--input_file_chain2', required=True)
    args.add_argument(f'--n_processes', default=8)

    return args.parse_args()


def get_tool_input_from_chain_dfs(chain_df1, chain_df2, translate_short_nt_to_cdr3_aa=True, v_dom_to_aa=True, scfv_dom_to_aa=True):
    chain_df1["chain"] = 1
    chain_df2["chain"] = 2

    chain_df1["row"] =  chain_df1["curated_receptor_id"].astype('Int64').astype(str) # "cur" + chain_df1["chain1_id_cur"].astype('Int64').astype(str) + "-calc" + chain_df1["chain1_id_calc"].astype('Int64').astype(str)
    chain_df2["row"] =  chain_df2["curated_receptor_id"].astype('Int64').astype(str) # "cur" + chain_df2["chain2_id_cur"].astype('Int64').astype(str) + "-calc" + chain_df2["chain2_id_calc"].astype('Int64').astype(str)

    rename_cols = {"v_gene": "vgene_curated", "d_gene": "dgene_curated", "j_gene": "jgene_curated",
                   "cdr1_seq": "cdr1_seq_curated", "cdr2_seq": "cdr2_seq_curated", "cdr3_seq": "cdr3_seq_curated",
                   "scientific_name": "species_latin"}

    chain_df1.rename(columns={**rename_cols,
                              **{"chain1_nt_seq": "nt",
                           "chain1_pro_seq": "aa",
                           "chain1_id_cur": "chain_id_cur",
                           "chain1_id_calc": "chain_id_calc"}}, inplace=True)

    chain_df2.rename(columns={**rename_cols,
                              **{"chain2_nt_seq": "nt",
                           "chain2_pro_seq": "aa",
                           "chain2_id_cur": "chain_id_cur",
                           "chain2_id_calc": "chain_id_calc"}}, inplace=True)


    assert len(chain_df1["row"].unique()) == len(chain_df1), "Duplicate chain ids"
    assert len(chain_df2["row"].unique()) == len(chain_df2), "Duplicate chain ids"

    df = pd.concat([chain_df1[util.TOOL_INPUT_COLUMNS], chain_df2[util.TOOL_INPUT_COLUMNS]])

    return format_tool_input_columns(df, translate_short_nt_to_cdr3_aa, v_dom_to_aa, scfv_dom_to_aa)


def omit_from_calculdated_data(row):
    with set_logging_context(template_row=row["row"], chain=row["chain"]):
        # Omit chains with incorrect V/D/J genes
        if ((pd.notna(row["v_gene"]) and row["v_gene"][3] != "V") or
                (pd.notna(row["d_gene"]) and row["d_gene"][3] != "D") or
                (pd.notna(row["j_gene"]) and row["j_gene"][3] != "J")):
            logging.error("Omitting from calculated data: V/D/J gene incorrect, re-curation may be needed")
            return True

        # Omit chains with illegal characters in CDRs
        for seq_col in ("cdr1_seq", "cdr2_seq", "cdr3_seq", "junction_seq"):
            if pd.notna(row[seq_col]):
                if not all(aa in util.AA_ALPHABET_STRICT for aa in set(row[seq_col])):
                    logging.error(f"Omitting from calculated data: illegal character in {seq_col}: {row[seq_col]}, re-curation may be needed")
                    return True

        # If the chain type changed, curated data updates are required
        if row["chain_type_calc"] != row["chain_type_cur"]:
            if pd.isna(row["chain_type_cur"]) and pd.isna(row["chain_type_calc"]):
                pass
            elif not (row["chain_type_cur"] == "light" and row["chain_type_calc"].endswith("_light")):
                logging.error(f"Omitting from calculated data: changed chain type (row: {row["row"]}, curated: {row["chain_type_cur"]}, calculated: {row["chain_type_calc"]}), re-curation may be needed")
                return True

    return False


def validate_tool_input(tool_input_full, output_folder):
    is_valid = tool_input_full.apply(validate_row_db, axis=1)
    invalid = tool_input_full[~is_valid]
    invalid.to_csv(output_folder / "tool_input_invalid.csv", index=False)

    return tool_input_full[is_valid].copy()


def write_db_recalculated_chains(tool_input, consolidate_results_df, output_folder):
    tool_input["chain_id_calc"] = tool_input["chain_id_calc"].astype("Int64")

    consolidate_results_df["row"] = consolidate_results_df["row"].astype("str")

    calculated_chains_df = pd.merge(consolidate_results_df,
                                    tool_input[["row", "chain", "chain_id_cur", "chain_id_calc", "chain_type"]],
                                    how="left", on=["row", "chain"], suffixes=("_calc", "_cur"))

    to_omit = calculated_chains_df.apply(omit_from_calculdated_data, axis=1)

    calculated_chains_df[to_omit].to_csv(output_folder / "recalculated_chains_failed.csv", index=False)

    calculated_chains_success_df = calculated_chains_df[~to_omit].copy()

    calculated_chains_success_df.rename(columns={"chain_id_calc": "chain_id", "chain_id_cur": "chain_id_cur_for_reference_only", "chain_type_calc": "chain_type"}, inplace=True)
    calculated_chains_success_df = calculated_chains_success_df[["chain_id_cur_for_reference_only",
                                                                 "chain_id", "species", "chain_type", "v_dom_seq",
                                                                 "v_gene", "d_gene", "j_gene",
                                                                 "cdr1_seq", "cdr2_seq", "cdr3_seq", "junction_seq"]]

    calculated_chains_success_df[calculated_chains_success_df["chain_id"].notna()].to_csv(output_folder / "recalculated_chains_success.csv", index=False)
    calculated_chains_success_df[calculated_chains_success_df["chain_id"].isna()].to_csv(output_folder / "recalculated_chains_success_no_identifier.csv", index=False)


def main(args):
    output_folder, tmp_folder, intermediate_folder = setup(args.output_folder)

    validate = True
    compute_results = True

    df1 = pd.read_csv(args.input_file_chain1)
    df2 = pd.read_csv(args.input_file_chain2)

    tool_input = get_tool_input_from_chain_dfs(df1, df2)

    if validate:
        tool_input = validate_tool_input(tool_input, output_folder)

    consolidate_results_df = run(intermediate_folder, tmp_folder, tool_input, compute_results)

    print("Writing recalculated chains files...")
    write_db_recalculated_chains(tool_input, consolidate_results_df, output_folder)
    print("...done")


if __name__ == '__main__':
    args = parse_arguments()
    main(args)