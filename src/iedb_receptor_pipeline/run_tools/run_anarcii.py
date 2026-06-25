import iedb_receptor_pipeline.utilities.util as util
from iedb_receptor_pipeline.utilities.logging_setup import set_logging_context, LOG_FILE_SEPARATOR
import pandas as pd
from anarcii import Anarcii
import logging


def get_sequence_from_domain(domain, region=None, imgt_range_strict=True, aa_strict=False):
    is_complete = True

    if region is not None:
        imgt_range = util.CDR_RANGES[region]

        if imgt_range_strict:
            if not all([idx in [item[0][0] for item in domain] for idx in imgt_range]):
                logging.info(f"the domain computed by ANARCII is incomplete and the full '{region}' field cannot be calculated.")
                is_complete = False

        domain = [item for item in domain if item[0][0] in imgt_range]

    seq = "".join([item[1] for item in domain if item[1] != "-"])
    seq = seq if is_complete and len(seq) > 0 else None

    if aa_strict and seq is not None:
        seq = seq if all(aa in util.AA_ALPHABET_STRICT for aa in seq) else None

    return seq

def get_anarcii_input(full_input):
    valid_chains = util.BCR_CHAINS.union(util.TCR_CHAINS).union({"ignar"})
    columns_of_interest = ["row", "chain", "species_latin", "aa", "chain_type", "receptor_type"]
    candidate_rows = full_input["chain_type"].isin(valid_chains) & full_input["aa"].notna() & (full_input["aa"].str.len() > 40)

    anarcii_input = full_input[candidate_rows][columns_of_interest].copy()
    anarcii_input["seq_type_anarci"] = anarcii_input["chain_type"].map(lambda x: "antibody" if x in util.BCR_CHAINS else "tcr")

    anarcii_input.loc[anarcii_input["species_latin"].isin(util.IEDB_SHARK_SPECIES), "seq_type_anarci"] = "vnar"
    anarcii_input.loc[anarcii_input["chain_type"].str.lower() == "ignar", "seq_type_anarci"] = "vnar"
    anarcii_input.loc[anarcii_input["species_latin"].isin(util.IEDB_CAMEL_SPECIES) & anarcii_input["receptor_type"].isin(["heavy", "heavyheavy"]), "seq_type_anarci"] = "vhh"

    return anarcii_input

def get_anarcii_annotations(full_input, n_processes):
    anarcii_input = get_anarcii_input(full_input)

    structured_results = []

    for (seq_type, chain_no, species), subset_input in anarcii_input.groupby(['seq_type_anarci', 'chain', 'species_latin']):
        model = Anarcii(seq_type=seq_type, mode="accuracy", ncpu=n_processes)
        results = model.number(dict(zip(subset_input["row"], subset_input["aa"])))
        structured_results.append(restructure_anarcii_output(results, chain_no, species))

    if len(structured_results) > 0:
        all_annotations = anarcii_input.merge(pd.concat(structured_results), on=["row", "chain"], suffixes=("_curated", "_anarcii"))
        all_annotations["chain_type"] = all_annotations["chain_type_anarcii"]
        all_annotations = all_annotations.drop(["chain_type_curated", "chain_type_anarcii"], axis=1)

        return all_annotations


def get_reason_junction_invalid(junction, v_domain, species):
    if not all(aa in util.AA_ALPHABET_STRICT for aa in junction):
        return "Junction contains illegal characters"

    if len(junction) < 6:
        if v_domain.endswith(junction):
            return "V domain too short"
        return f"Junction is too short ({len(junction)})"

    if junction[0] != "C":
        return "Conserved starting C not found"

    if species in ("Homo sapiens", "Mus musculus"):
        legal_trailing_aas = ("F", "W", "C") if species == "Homo sapiens" else ("F", "W", "L", "C")

        if junction[-1] not in legal_trailing_aas:
            if v_domain.endswith(junction):
                return "V domain too short"

            return f"Conserved trailing aa ({legal_trailing_aas}) not found"

    return None


def restructure_anarcii_output(results, chain_no, species):
    df = pd.DataFrame(columns=["row", "chain", "chain_type", "v_dom_seq",
                               "cdr1_seq_calculated", "cdr1_seq_complete",
                               "cdr2_seq_calculated", "cdr2_seq_complete",
                               "cdr3_seq_calculated", "cdr3_seq_complete",
                               "junction_seq_calculated", "junction_seq_complete",
                               "vgene_calculated", "jgene_calculated"])

    for row, row_result in results.items():
        with set_logging_context(template_row=row, chain=chain_no):

            if row_result["error"] is not None:
                logging.error(f"ANARCII error: {row_result["error"]}")

            if row_result["numbering"] is None:
                continue

            v_domain = get_sequence_from_domain(row_result["numbering"], region=None)
            cdr1_seq_calculated = get_sequence_from_domain(row_result["numbering"], region="cdr1", aa_strict=True)
            cdr2_seq_calculated = get_sequence_from_domain(row_result["numbering"], region="cdr2", aa_strict=True)
            junction_seq_calculated = get_sequence_from_domain(row_result["numbering"], region="junction", aa_strict=True)

            cdr3_seq_calculated = None

            if junction_seq_calculated is not None:
                invalid_reason = get_reason_junction_invalid(junction_seq_calculated, v_domain, species)

                if invalid_reason is not None:
                    if invalid_reason == "V domain too short":
                        logging.info(f"Failed to compute a valid ANARCII CDR3/Junction from {species} {junction_seq_calculated}: {invalid_reason}{LOG_FILE_SEPARATOR}full V domain: {v_domain}")
                    else:
                        logging.warning(f"Failed to compute a valid ANARCII CDR3/Junction from {species} {junction_seq_calculated}: {invalid_reason}{LOG_FILE_SEPARATOR}full V domain: {v_domain}")

                    junction_seq_calculated = None

            if junction_seq_calculated is not None:
                cdr3_seq_calculated = junction_seq_calculated[1:-1]

            df.loc[len(df)] = {"row": row,
                               "chain": chain_no,
                               "chain_type": util.CHAIN_ANARCI_TO_IEDB_SINGLE[row_result["chain_type"]],
                               "cdr1_seq_calculated": cdr1_seq_calculated,
                               "cdr2_seq_calculated": cdr2_seq_calculated,
                               "cdr3_seq_calculated": cdr3_seq_calculated,
                               "junction_seq_calculated": junction_seq_calculated,
                               "v_dom_seq": v_domain}

    return df
