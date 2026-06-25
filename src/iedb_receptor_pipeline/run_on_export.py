from argparse import ArgumentParser


from parse_input import get_tool_input
from iedb_receptor_pipeline.run_tools.run_tools import run, setup

import numpy as np
from iedb_receptor_pipeline.utilities import util
import pandas as pd

bcr_tcr = "tcr"


def parse_arguments():
    args = ArgumentParser()

    args.add_argument(f'--input_file', required=True)
    args.add_argument(f'--output_folder',  default=f"./iedb_receptors_output")

    return args.parse_args()


def overwrite_calc_data(full_input_orig, consolidate_results_df):
    full_input_orig.drop(columns=[('Chain 1', 'CDR1 Start Curated'), ('Chain 1', 'CDR1 End Curated'),
                                  ('Chain 1', 'CDR2 Start Curated'), ('Chain 1', 'CDR2 End Curated'),
                                  ('Chain 1', 'CDR2 Start Calculated'), ('Chain 1', 'CDR2 End Calculated'),
                                  ('Chain 1', 'CDR1 Start Calculated'), ('Chain 1', 'CDR1 End Calculated'),
                                  ('Chain 1', 'CDR3 Start Curated'), ('Chain 1', 'CDR3 End Curated'),
                                  ('Chain 1', 'CDR3 Start Calculated'), ('Chain 1', 'CDR3 End Calculated'),
                                  ('Chain 2', 'CDR3 Start Calculated'), ('Chain 2', 'CDR3 End Calculated'),
                                  ('Chain 2', 'CDR3 Start Curated'), ('Chain 2', 'CDR3 End Curated'),
                                  ('Chain 2', 'CDR1 Start Curated'), ('Chain 2', 'CDR1 End Curated'),
                                  ('Chain 2', 'CDR1 Start Calculated'), ('Chain 2', 'CDR1 End Calculated'),
                                  ('Chain 2', 'CDR2 Start Curated'), ('Chain 2', 'CDR2 End Curated'),
                                  ('Chain 2', 'CDR2 Start Calculated'), ('Chain 2', 'CDR2 End Calculated')],
                         inplace=True)

    full_input_orig[('Chain 1', 'Junction Calculated')] = np.nan
    full_input_orig[('Chain 2', 'Junction Calculated')] = np.nan
    full_input_orig[('Chain 1', 'V Domain Calculated')] = np.nan
    full_input_orig[('Chain 2', 'V Domain Calculated')] = np.nan

    for i in full_input_orig.index:
        if i % 100 == 0:
            print(i, "/", len(full_input_orig))
        for chain in (1, 2):
            calc_row = consolidate_results_df[(consolidate_results_df["row"] == i) & (consolidate_results_df["chain"] == chain)]
            assert len(calc_row) <= 1, f"len = {len(calc_row)}"

            if len(calc_row) == 1:
                full_input_orig.loc[i, (f'Chain {chain}', 'Type')] = calc_row["chain_type"].tolist()[0]
                full_input_orig.loc[i, (f'Chain {chain}', 'Calculated V Gene')] = calc_row["v_gene"].tolist()[0]
                full_input_orig.loc[i, (f'Chain {chain}', 'Calculated D Gene')] = calc_row["d_gene"].tolist()[0]
                full_input_orig.loc[i, (f'Chain {chain}', 'Calculated J Gene')] = calc_row["j_gene"].tolist()[0]
                full_input_orig.loc[i, (f'Chain {chain}', 'CDR1 Calculated')] = calc_row["cdr1_seq"].tolist()[0]
                full_input_orig.loc[i, (f'Chain {chain}', 'CDR2 Calculated')] = calc_row["cdr2_seq"].tolist()[0]
                full_input_orig.loc[i, (f'Chain {chain}', 'CDR3 Calculated')] = calc_row["cdr3_seq"].tolist()[0]
                full_input_orig.loc[i, (f'Chain {chain}', 'Junction Calculated')] = calc_row["junction_seq"].tolist()[0]
                full_input_orig.loc[i, (f'Chain {chain}', 'V Domain Calculated')] = calc_row["v_dom_seq"].tolist()[0]


    return full_input_orig.copy()


def main(args):
    output_folder, tmp_folder, intermediate_folder = setup(args.output_folder)

    full_input_orig = pd.read_csv(args.input_file, header=[0, 1], sep=",", low_memory=False)
    full_input = full_input_orig[list(util.DOUBLE_TO_SINGLE_HEADER.keys())].copy()
    full_input.columns = [util.DOUBLE_TO_SINGLE_HEADER[col] for col in full_input.columns.values]

    tool_input = get_tool_input(full_input, retrieve_accessions=False)

    compute_results = True

    consolidate_results_df = run(output_folder, tmp_folder, tool_input, compute_results)

    overwrite_calc = True

    if overwrite_calc:
        full_output = overwrite_calc_data(full_input_orig, consolidate_results_df)
        full_output.to_csv(output_folder / f"{bcr_tcr}_full_v4.csv", index=False)




if __name__ == '__main__':
    args = parse_arguments()
    main(args)
