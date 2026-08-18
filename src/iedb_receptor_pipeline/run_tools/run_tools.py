import iedb_receptor_pipeline.utilities.util as util
from iedb_receptor_pipeline.run_tools.run_anarcii import get_anarcii_annotations
from iedb_receptor_pipeline.run_tools.run_igblast import run_igblast
from iedb_receptor_pipeline.run_tools.run_tidytcells import fix_with_tidytcells
from iedb_receptor_pipeline.consolidate_results import consolidate_results

from pathlib import Path
import logging

import pandas as pd




def setup(output_folder):
    output_folder = Path(output_folder)
    intermediate_folder = output_folder / "intermediate_results"
    tmp_folder = output_folder / "tmp"

    util.build_path(output_folder, delete_if_exists=False)
    util.build_path(intermediate_folder, delete_if_exists=True)
    util.build_path(tmp_folder, delete_if_exists=True)

    logging.basicConfig(filename=output_folder / "log.txt", level=logging.WARNING,
                        format="%(asctime)s %(levelname)s:%(message)s")


    return output_folder, tmp_folder, intermediate_folder



def run(output_folder, tmp_folder, tool_input, compute_results, n_processes=8):
    tool_input.to_csv(output_folder / "tool_input.csv")

    if compute_results:
        print("running tt...")
        tt_results = fix_with_tidytcells(tool_input)
        tt_results.to_csv(output_folder / "tt_results.csv")
        print("tt done")

        print("running anarcii...")
        anarcii_results = get_anarcii_annotations(tool_input, n_processes)

        if anarcii_results is not None:
            anarcii_results.to_csv(output_folder / "anarcii_results.csv")
        mssg = ", no results computed" if anarcii_results is None else ""
        print("anarcii done" + mssg)

        print("running igblast...")
        igblast_nt_results, igblast_aa_results = run_igblast(tool_input, tmp_folder, n_processes)
        if igblast_nt_results is not None:
            igblast_nt_results.to_csv(output_folder / "igblast_nt_results.csv")

        if igblast_aa_results is not None:
            igblast_aa_results.to_csv(output_folder / "igblast_aa_results.csv")
        mssg = ", no results computed" if (igblast_nt_results is None and igblast_aa_results is None) else ""
        print("igblast done" + mssg)


    tool_input = pd.read_csv(output_folder / "tool_input.csv")
    tt_results = util.safe_read_csv(output_folder / "tt_results.csv")
    anarcii_results = util.safe_read_csv(output_folder / "anarcii_results.csv")
    igblast_nt_results = util.safe_read_csv(output_folder / "igblast_nt_results.csv")
    igblast_aa_results = util.safe_read_csv(output_folder / "igblast_aa_results.csv")


    consolidate_results_df = consolidate_results(tool_input, anarcii_results, igblast_nt_results, igblast_aa_results, tt_results)
    consolidate_results_df.to_csv(output_folder / "consolidate_results.csv")

    return consolidate_results_df

