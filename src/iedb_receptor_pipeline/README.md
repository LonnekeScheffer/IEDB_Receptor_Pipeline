# IEDB receptor data processing pipeline overview

The processing pipeline has the following steps:

1. The pipeline can be run using various entrypoints:
    - [run_on_export.py](run_on_export.py) runs on the full IEDB TCR/BCR export files
    - [run_on_curation_template.py](run_on_curation_template.py) runs on the IEDB curation template for new data
    - [run_on_db.py](run_on_db.py) runs on direct database exports


2. All entrypoints use sanity checks in [parse_input.py](parse_input.py) for data validation. 


3. Next, [run_tools.py](run_tools/run_tools.py) will execute [run_anarcii.py](run_tools/run_anarcii.py), [run_igblast.py](run_tools/run_igblast.py) and [run_tidytcells.py](run_tools/run_tidytcells.py) to produce calculated data outputs. 


4. Finally, [consolidate_results.py](consolidate_results.py) combines the data from the different tools, and performs the final sanity checks. 

