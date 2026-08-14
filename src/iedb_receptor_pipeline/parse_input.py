from iedb_receptor_pipeline.utilities import util
from iedb_receptor_pipeline.utilities.logging_setup import set_logging_context, LOG_FILE_SEPARATOR
import logging
import pandas as pd

def get_illegal_chars(seq, alphabet):
    return "'" + "', '".join(list(set(seq).difference(alphabet))) + "'"


def warn_length(seq, column, MIN_LENGTH, MAX_LENGTH_INFO, MAX_LENGTH_WARN, seqtype, is_scfv=False, extra_info="", skip_max_length_warning=True):
    scfv_info = "for scFv sequences" if is_scfv else ""

    if pd.notna(seq):
        if len(seq) <= MIN_LENGTH:
            if seqtype == "nucleotide" and len(seq) <= 40*3:
                logging.info(f"Length of {seqtype} sequence in {column} is {len(seq)}, an attempt will be made to translate this to the CDR3. ")
            else:
                logging.warning(f"Length of {seqtype} sequence in {column} is {len(seq)}, this is shorter than expected {scfv_info} (<= {MIN_LENGTH} {seqtype}s){LOG_FILE_SEPARATOR} {seq} {extra_info}")

        # too long sequences are very common (C region included)
        if len(seq) >= MAX_LENGTH_INFO and skip_max_length_warning == False:
            if len(seq) >= MAX_LENGTH_WARN:
                logging.warning(f"Length of {seqtype} sequence in {column} is {len(seq)}, this is longer than expected {scfv_info} (>= {MAX_LENGTH_WARN} {seqtype}s){LOG_FILE_SEPARATOR} {seq} {extra_info}")
            else:
                logging.info(f"Length of {seqtype} sequence in {column} is {len(seq)}, this is longer than expected {scfv_info} (>= {MAX_LENGTH_INFO} {seqtype}s){LOG_FILE_SEPARATOR} {seq} {extra_info}")


def validate_row_chain(row, receptor_type_col, aa_col, v_dom_col, cdr1_col, cdr2_col, cdr3_col, nt_col, v_gene_col, d_gene_col, j_gene_col, chain_type_col, species_col):
    valid = True

    FULL_AA_MIN_LENGTH = 85
    FULL_AA_MAX_LENGTH = 170

    SCFV_FULL_AA_MIN_LENGTH = FULL_AA_MIN_LENGTH * 2
    SCFV_FULL_AA_MAX_LENGTH = FULL_AA_MAX_LENGTH * 2 + 20

    if pd.isna(receptor_type_col):
        logging.error("Receptor type is missing")
        valid = False

    if all([pd.isna(row[field]) for field in [aa_col, v_dom_col, cdr1_col, cdr2_col, cdr3_col, nt_col, v_gene_col, d_gene_col, j_gene_col]]):
        # all data is NA, chain and species should be NA as well
        if pd.notna(row[chain_type_col]) or pd.notna(row[species_col]):
            logging.warning("This chain has no data, please remove chain type and species")
    else:
        # Some data is available, chain and species type must be present
        if pd.isna(row[chain_type_col]):
            logging.error("Chain type is missing")
            valid = False

        if pd.isna(row[species_col]):
            logging.error("Species is missing")
            valid = False

    try:
        resolve_species(row[species_col])
    except AssertionError:
        logging.error(f"Unsupported {species_col} ontology: {row[species_col]}")
        valid = False

    for aa_seq_expected_column in [cdr1_col, cdr2_col, cdr3_col]:
        if aa_seq_expected_column is not None:
            if pd.notna(row[aa_seq_expected_column]) and not util.is_valid_alphabet(row[aa_seq_expected_column].strip(), util.AA_ALPHABET_STRICT):
                logging.error(f"Illegal characters in {aa_seq_expected_column}: {get_illegal_chars(row[aa_seq_expected_column], util.AA_ALPHABET_STRICT)}{LOG_FILE_SEPARATOR}full_sequence: {row[aa_seq_expected_column]}")
                valid = False

    for aa_seq_expected_column in [aa_col, v_dom_col]:
        if aa_seq_expected_column is not None:
            if pd.notna(row[aa_seq_expected_column]) and not util.is_valid_alphabet(row[aa_seq_expected_column].strip(), util.AA_ALPHABET_AMBIGUOUS):
                logging.error(f"Illegal characters in {aa_seq_expected_column}: {get_illegal_chars(row[aa_seq_expected_column], util.AA_ALPHABET_AMBIGUOUS)}{LOG_FILE_SEPARATOR}full_sequence: {row[aa_seq_expected_column]}")
                valid = False

    warn_length(row[cdr1_col], cdr1_col, 4, 15, 15, "amino acid", extra_info=f"(V gene: {row[v_gene_col]}, J gene: {row[j_gene_col]}, full seq available: {pd.notna(row[aa_col])})")
    warn_length(row[cdr2_col], cdr2_col, 3, 15, 15, "amino acid", extra_info=f"(V gene: {row[v_gene_col]}, J gene: {row[j_gene_col]}, full seq available: {pd.notna(row[aa_col])}")
    warn_length(row[cdr3_col], cdr3_col, 4, 35, 35, "amino acid", extra_info=f"(V gene: {row[v_gene_col]}, J gene: {row[j_gene_col]}, full seq available: {pd.notna(row[aa_col])}")

    if v_dom_col is not None:
        warn_length(row[v_dom_col], v_dom_col, FULL_AA_MIN_LENGTH, FULL_AA_MAX_LENGTH, FULL_AA_MAX_LENGTH, "amino acid")

    if aa_col is not None and pd.notna(row[aa_col]):
        if str(row[receptor_type_col]).lower() in ("scfv", "tscfv") and util.is_scfv(row[aa_col]):
            warn_length(row[aa_col], aa_col, FULL_AA_MIN_LENGTH, SCFV_FULL_AA_MAX_LENGTH, SCFV_FULL_AA_MAX_LENGTH*2, f"amino acid", is_scfv=True)
        else:
            if (str(row[receptor_type_col]).lower() not in ("scfv", "tscfv")) and (len(row[aa_col]) > FULL_AA_MIN_LENGTH) and util.is_scfv(row[aa_col]):
                logging.warning(f"Long sequence with possible SG-linker detected, this receptor might be a (T)scFv but receptor type was {str(row[receptor_type_col]).lower()}")

            warn_length(row[aa_col], aa_col, FULL_AA_MIN_LENGTH, FULL_AA_MAX_LENGTH, FULL_AA_MAX_LENGTH*2, "amino acid")

    if v_dom_col is not None and pd.notna(row[aa_col]) and pd.notna(row[v_dom_col]): # todo same for cdr3/cdr2/etc?
        if row[v_dom_col][1:] not in row[aa_col]:
            logging.error(f"Expected V domain to be a subsequence of the full aa sequence{LOG_FILE_SEPARATOR}(V domain: {row[v_dom_col]}; full sequence: {row[aa_col]})")
            #valid = False

    if pd.notna(row[nt_col]):
        if not util.is_valid_alphabet(row[nt_col].lower().strip(), util.NT_ALPHABET):
            if not util.is_valid_alphabet(row[nt_col].lower().strip(), util.NT_ALPHABET_AMBIGUOUS):
                valid = False
                logging.error(f"Illegal characters in {nt_col}: {get_illegal_chars(row[nt_col].lower().strip(), util.NT_ALPHABET)}{LOG_FILE_SEPARATOR}full sequence: {row[nt_col].lower().strip()}")
            else:
                logging.info(f"Illegal characters in {nt_col}: {get_illegal_chars(row[nt_col].lower().strip(), util.NT_ALPHABET)}{LOG_FILE_SEPARATOR}full sequence: {row[nt_col].lower().strip()}")

        if len(row[nt_col]) % 3 != 0:
            logging.info(f"Length of nt sequence in {nt_col} is not a multiple of 3.")

        if str(row[receptor_type_col]).lower() in ("scfv", "tscfv"):
            warn_length(row[nt_col], nt_col, SCFV_FULL_AA_MIN_LENGTH * 3, SCFV_FULL_AA_MAX_LENGTH * 3, SCFV_FULL_AA_MAX_LENGTH * 6, "nucleotide", is_scfv=True)
        else:
            warn_length(row[nt_col], nt_col, FULL_AA_MIN_LENGTH * 3, FULL_AA_MAX_LENGTH * 3, FULL_AA_MAX_LENGTH * 6, "nucleotide")

    if pd.notna(row[chain_type_col]):
        if row[chain_type_col] in util.VJ_CHAINS:
            if pd.notna(row[d_gene_col]):
                logging.error(f"Chain type {row[chain_type_col]} is not expected to have D gene.")
                valid = False

    if str(row[receptor_type_col]).lower() not in util.RECEPTOR_TYPE_TO_VALID_CHAINS.keys():
        logging.error(f"Unrecognized receptor type '{str(row[receptor_type_col]).lower()}' (expected one of {list(util.RECEPTOR_TYPE_TO_VALID_CHAINS.keys())}.")
        valid = False
    else:
        if pd.notna(row[chain_type_col]):
            if str(row[chain_type_col]).lower() not in util.RECEPTOR_TYPE_TO_VALID_CHAINS[str(row[receptor_type_col]).lower()]:
                logging.error(f"Unrecognized chain type {row[chain_type_col]} for receptor type {row[receptor_type_col]} (expected one of {util.RECEPTOR_TYPE_TO_VALID_CHAINS[str(row[receptor_type_col]).lower()]}.")
                valid = False

    return valid

def get_name_from_ontology(ontology_source, ontology_number):
    if ontology_source.upper() == "NCBITAXON":
        return util.NCBITAXON_TO_LATIN[ontology_number]
    elif ontology_source.upper() == "ONTIE":
        return util.ONTIE_TO_LATIN[ontology_number]

    if ontology_number in util.NCBITAXON_TO_LATIN:
        return util.NCBITAXON_TO_LATIN[ontology_number]
    elif ontology_number in util.ONTIE_TO_LATIN:
        return util.ONTIE_TO_LATIN[ontology_number]
    elif ontology_number in util.ORGANISM_ID_TO_LATIN:
        return util.ORGANISM_ID_TO_LATIN[ontology_number]

    assert False, f"Cannot resolve ontology from {ontology_source} {ontology_number}"
    return "unknown"


def resolve_species(species_input):
    if pd.isna(species_input):
        return None

    if species_input in util.NCBITAXON_TO_LATIN.values():
        return species_input

    if type(species_input) in (int, float):
        return get_name_from_ontology("unknown", str(int(species_input)))

    if type(species_input) == str:
        if species_input != species_input.strip():
            logging.warning(f"Whitespace characters were removed from species input, please correct to remove this warning: '{species_input}'")
            species_input = species_input.strip()

        if species_input.isdigit():
            return get_name_from_ontology("unknown", species_input)

        if "/" in species_input:
            ontology_str = species_input.rsplit("/", 1)[1]

            if "_"  in ontology_str:
                ontology_source, ontology_number = ontology_str.rsplit("_")
            elif ":" in ontology_str:
                ontology_source, ontology_number = ontology_str.rsplit(":")
            else:
                assert False, f"Cannot resolve ontology from {species_input}"
        else:
            assert False, f"Cannot resolve ontology from {species_input}"

        return get_name_from_ontology(ontology_source, ontology_number)

    assert False, f"Cannot resolve ontology from {species_input}"

def get_chain_columns(curator_df, chain):
    chain_df = pd.DataFrame(columns=["chain", "species_latin", "aa", "nt"], index=curator_df.index)

    chain_df["chain"] = chain
    chain_df["curated_receptor_id"] = curator_df[f"curated_receptor_id"]
    chain_df["chain_id_cur"] = curator_df[f"curated_chain{chain}_id"]
    chain_df["chain_id_calc"] = curator_df[f"calculated_chain{chain}_id"]
    chain_df["receptor_type"] = curator_df[f"receptor_type"]
    chain_df["chain_type"] = curator_df[f"chain{chain}_type"]
    chain_df["species_latin"] = curator_df[f"chain{chain}_species"].apply(resolve_species)

    chain_df["organism_name"] = chain_df["species_latin"] # todo distinguish between species and organism? -> for mouse strain-specific gene annotation

    chain_df["aa"] = curator_df[f"chain{chain}_pro_seq"]
    chain_df["nt"] = curator_df[f"chain{chain}_nucleotide"]

    chain_df["v_dom_seq"] = curator_df[f"chain{chain}_vdomain_curated"]

    chain_df["vgene_curated"] = curator_df[f"chain{chain}_v_gene_curated"]
    chain_df["dgene_curated"] = curator_df[f"chain{chain}_d_gene_curated"]
    chain_df["jgene_curated"] = curator_df[f"chain{chain}_j_gene_curated"]

    chain_df["cdr1_seq_curated"] = curator_df[f"chain{chain}_cdr1_seq_curated"]
    chain_df["cdr2_seq_curated"] = curator_df[f"chain{chain}_cdr2_seq_curated"]
    chain_df["cdr3_seq_curated"] = curator_df[f"chain{chain}_cdr3_seq_curated"]

    chain_df = chain_df.reset_index().rename(columns={"index": "row"})
    chain_df["row"] = chain_df["row"].astype(str)

    return chain_df

def validate_ids(row):
    for id_col in ["ASSAY_ID", "REFERENCE_ID"]:
        if id_col in row:
            if pd.isna(row[id_col]) or row[id_col] == "":
                logging.error(f"{id_col} not found")
            else:
                for id in row[id_col].split(", "):
                    if not id.strip().isnumeric():
                        logging.error(f"{id_col} expected numeric identifier: {id}")

                if "," in row[id_col] and (len(row[id_col].split(", ")) != len(set(row[id_col].split(", ")))):
                    logging.error(f"{id_col} expected unique identifiers, found duplicates: {row[id_col]}")


def validate_row_curation_template(row):
    all_chains_valid = []

    for chain in (1, 2):
        with set_logging_context(template_row=row.name, chain=chain):
            validate_ids(row)

            is_valid = validate_row_chain(row=row,
                                          receptor_type_col="receptor_type",
                                          aa_col=f"chain{chain}_pro_seq",
                                          v_dom_col=f"chain{chain}_vdomain_curated",
                                          cdr1_col=f"chain{chain}_cdr1_seq_curated",
                                          cdr2_col=f"chain{chain}_cdr2_seq_curated",
                                          cdr3_col=f"chain{chain}_cdr3_seq_curated",
                                          nt_col=f"chain{chain}_nucleotide",
                                          chain_type_col=f"chain{chain}_type",
                                          v_gene_col=f"chain{chain}_v_gene_curated",
                                          d_gene_col=f"chain{chain}_d_gene_curated",
                                          j_gene_col=f"chain{chain}_j_gene_curated",
                                          species_col=f"chain{chain}_species")

            all_chains_valid.append(is_valid)

    return all(all_chains_valid)


def validate_row_db(row):
    with set_logging_context(template_row=row["row"], chain=row["chain"]):
        # if row["species_latin"] != row["organism_name"]:
        #     if row["organism_name"] in {''}:
        #         pass
        # todo deal with scfv

        return validate_row_chain(row=row,
                                  receptor_type_col="receptor_type",
                                  aa_col="aa",
                                  v_dom_col="v_dom_seq",
                                  cdr1_col="cdr1_seq_curated",
                                  cdr2_col="cdr2_seq_curated",
                                  cdr3_col="cdr3_seq_curated",
                                  nt_col="nt",
                                  v_gene_col="vgene_curated",
                                  d_gene_col="dgene_curated",
                                  j_gene_col="jgene_curated",
                                  chain_type_col="chain_type",
                                  species_col="species_latin")

def safe_translate(nt):
    from Bio.Seq import Seq

    try:
        nt = nt.lower().replace(" ", "")
        if len(nt) % 3 == 0:
            return str(Seq(nt).translate())
    except Exception as e:
        pass

def format_tool_input_columns(df, translate_short_nt_to_cdr3_aa, v_dom_to_aa, scfv_dom_to_aa):
    assert set(df.columns) == set(util.TOOL_INPUT_COLUMNS), f"df does not have expected columns: {df.columns}"

    df["row"] = df["row"].astype(str)

    # store the original input 'pro_seq', some modifications may be made to 'aa'
    df["pro_seq_orig"] = df["aa"]

    # if AA full seq is unavailable, but V_domain is available -> use V_domain instead
    if v_dom_to_aa:
        mask = df["aa"].isna() & df["v_dom_seq"].notna()
        df.loc[mask, "aa"] = df.loc[mask, "v_dom_seq"]

    #remove spaces from nt sequence
    df.loc[df["nt"].notna(), "nt"] = df.loc[df["nt"].notna(), "nt"].apply(lambda x: x.replace(" ", ""))

    # if AA full seq is unavailable, but NT is available -> translate (only used for CDR3)
    if translate_short_nt_to_cdr3_aa:
        if any(df["nt"].notna()):
            mask = df["nt"].notna() & (df["nt"].str.len() < 40) & df["cdr3_seq_curated"].isna()
            df.loc[mask, "cdr3_seq_curated"] = df.loc[mask, "nt"].apply(lambda x: safe_translate(x))

    scfvs = df[df["receptor_type"].str.lower().isin(["scfv", "tscfv"])]

    if len(scfvs) > 0:
        dual_chain_scfvs = scfvs[scfvs["row"].map(scfvs["row"].value_counts()) == 2]

        if len(dual_chain_scfvs) > 0:
            dual_chain_scfvs["aa"].apply(util.split_scfv)

            dual_chain_scfvs.loc[:, ["chain1_scfv_aa", "chain2_scfv_aa"]] = pd.DataFrame(
                dual_chain_scfvs["aa"].apply(util.split_scfv).tolist(),
                index=dual_chain_scfvs.index,
                columns=["chain1_scfv_aa", "chain2_scfv_aa"]
            )

            for chain in (1, 2):
                scfv_new_aa = dual_chain_scfvs[(dual_chain_scfvs["chain"] == chain)][["row", "chain", f"chain{chain}_scfv_aa"]].dropna(how="any")

                df = df.merge(scfv_new_aa, left_on=['row', 'chain'], right_on=['row', 'chain'], how='left')
                df["aa"] = df[f"chain{chain}_scfv_aa"].combine_first(df['aa'])

            df.drop(columns=["chain1_scfv_aa", "chain2_scfv_aa"], inplace=True)

        # if the AA full seq of scFvs is identical for chain 1 and chain 2, and V domains exist, overwrite with V domain
        if scfv_dom_to_aa:
            scfvs = df[df["receptor_type"].str.lower().isin(["scfv", "tscfv"]) & df["v_dom_seq"].notna()]
            if len(scfvs) > 0:
                pro_seq_comparison = scfvs.pivot(index="curated_receptor_id", columns="chain", values="pro_seq_orig").dropna().reset_index()
                scfvs_to_split = pro_seq_comparison[pro_seq_comparison[1] == pro_seq_comparison[2]]["curated_receptor_id"]
                mask = df["curated_receptor_id"].isin(scfvs_to_split)
                df.loc[mask, "aa"] = df.loc[mask, "v_dom_seq"]

    return df


def get_tool_input(curation_template, translate_short_nt_to_cdr3_aa=True, v_dom_to_aa=True, scfv_dom_to_aa=True):
    is_valid = curation_template.apply(validate_row_curation_template, axis=1)
    validated_input = curation_template[is_valid]

    chain1_df = get_chain_columns(validated_input, chain=1)
    chain2_df = get_chain_columns(validated_input, chain=2)

    df = format_tool_input_columns(pd.concat([chain1_df, chain2_df]), translate_short_nt_to_cdr3_aa, v_dom_to_aa, scfv_dom_to_aa)

    return df
