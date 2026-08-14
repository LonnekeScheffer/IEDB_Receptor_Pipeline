import tidytcells as tt
import iedb_receptor_pipeline.utilities.util as util
from iedb_receptor_pipeline.utilities.logging_setup import set_logging_context
import logging
import pandas as pd


def get_tr_ig_package(chain_type):
    if chain_type in util.BCR_CHAINS:
        return tt.ig
    elif chain_type in util.TCR_CHAINS:
        return tt.tr


def test_without_species(gene):
    if gene[0:3] in util.CHAIN_LOCUS_TO_IEDB.keys():
        return True
    logging.error(f"Testing gene '{gene}' in the absence of species. Not a valid locus: {gene[0:3]}")
    return False


def standardize_gene_single(module, gene, species):
    if species.replace(" ", "").lower() in ("homosapiens", "musmusculus"):
        result = module.standardize(symbol=gene, species=species, log_failures=True, allow_subgroup=True, enforce_functional=False)

        if result.is_standardized:
            return result.symbol

        logging.error(f"Failed to standardize {species} gene {gene}: {result.error}")

    else:
        logging.warning(f"Species '{species}' does not have an available gene set for gene '{gene}', testing in the absence of species..")
        if test_without_species(gene):
            return gene


def safe_standardize_gene(module, curated_gene_str, species):
    if pd.notna(curated_gene_str):
        if util.GENE_SEP in curated_gene_str:
            genes = curated_gene_str.split(util.GENE_SEP)
        elif " or " in curated_gene_str:
            genes = curated_gene_str.split(" or ")
        elif "," in curated_gene_str:
            genes = curated_gene_str.split(",")
        else:
            genes = [curated_gene_str]

        std_genes = []

        for gene in genes:
            result = standardize_gene_single(module, gene.strip(), species)

            if result is not None and result not in std_genes:
                std_genes.append(result)

        return util.GENE_SEP.join(std_genes)


def get_input_cdr3(row):
    input_cdr3 = None
    is_curated_cdr3 = None

    if pd.notna(row["cdr3_seq_curated"]):
        input_cdr3 = row["cdr3_seq_curated"]
        is_curated_cdr3 = True

    # If 'aa' contains a short sequence, presume this is CDR3
    elif pd.notna(row["cdr3_seq_curated"]) and (len(row["aa"]) < 40):
        input_cdr3 = row["aa"]
        is_curated_cdr3 = False

    return input_cdr3, is_curated_cdr3


def safe_return_junction_with_jgenes(junctions_with_jgenes, input_cdr3):
    junction_calc, junction_jgenes = list(junctions_with_jgenes.items())[0]

    if junction_jgenes == "" or junction_jgenes is None:
        junction_jgenes = None
    else:
        jgene_set_no_alleles = sorted(set({gene.rsplit("*")[0] for gene in junction_jgenes.split(", ")}))

        if len(jgene_set_no_alleles) > 1:
            logging.warning(
                f"J gene standardization of '{input_cdr3}' to '{junction_calc}' was based on multiple J gene matches ({junction_jgenes}), "
                f"The J annotations will be omitted.")
            junction_jgenes = None

        else:
            junction_jgenes = ", ".join(jgene_set_no_alleles)

    return junction_calc, junction_jgenes

def tt_standardize_junction(row, vgene_calc, jgene_calc, species, input_cdr3, is_curated_cdr3):
    junctions_with_jgenes = {}

    v_symbols = vgene_calc.split(util.GENE_SEP) if type(vgene_calc) == str and vgene_calc != "" else [None]
    j_symbols = jgene_calc.split(util.GENE_SEP) if type(jgene_calc) == str and jgene_calc != "" else [None]

    for v_symbol in v_symbols:
        for j_symbol in j_symbols:
            junction_calc, junction_jgene = tt_standardize_junction_single_vj(row, v_symbol, j_symbol, species, input_cdr3, is_curated_cdr3)

            if junction_calc is None and j_symbol is not None:
                junction_calc, junction_jgene = tt_standardize_junction_single_vj(row, v_symbol, None, species, input_cdr3, is_curated_cdr3)

            if junction_calc is not None and abs(len(junction_calc) - len(input_cdr3)) >= 6:
                logging.warning(f">=6 amino acids length difference between original CDR3 {input_cdr3} (v={v_symbol}, j={j_symbol}) and tidytcells-corrected junction {junction_calc}")

            if junction_calc is not None:
                junction_jgene = "" if junction_jgene is None else junction_jgene

                if junction_calc not in junctions_with_jgenes:
                    junctions_with_jgenes[junction_calc] = junction_jgene
                else:
                    junctions_with_jgenes[junction_calc] += util.GENE_SEP + junction_jgene

    if len(junctions_with_jgenes) > 1:
        logging.warning(f"Found multiple possible junction corrections based on different V/J genes: {junctions_with_jgenes.keys()}")

    elif len(junctions_with_jgenes) == 1:
        return safe_return_junction_with_jgenes(junctions_with_jgenes, input_cdr3)

    return None, None




def get_locus(chain_type):
    if chain_type == "heavy":
        return "IGH"
    if chain_type == "light":
        return "IG"
    elif chain_type== "kappa_light":
        return "IGK"
    elif chain_type== "lambda_light":
        return "IGL"
    elif chain_type== "alpha":
        return "TRA"
    elif chain_type== "beta":
        return "TRB"
    elif chain_type== "gamma":
        return "TRG"
    elif chain_type== "delta":
        return "TRD"

    assert False, "unknown chain type"

def tt_standardize_junction_single_vj(row, v_symbol, j_symbol, species, input_cdr3, is_curated_cdr3):
    if pd.isna(input_cdr3):
        return None, None

    # First: try only with functional genes
    std_junction = tt.junction.standardize(seq=input_cdr3,
                                           locus=get_locus(row['chain_type']),
                                           j_symbol=j_symbol,
                                           v_symbol=v_symbol,
                                           species=species,
                                           allow_c_correction=True if v_symbol is not None else False,
                                           allow_fw_correction=True if j_symbol is not None else False,
                                           max_v_reconstruction=3 if v_symbol is not None and is_curated_cdr3 else 1,
                                           max_j_reconstruction=3 if j_symbol is not None and is_curated_cdr3 else 1,
                                           enforce_functional_v=True,
                                           enforce_functional_j=True,
                                           log_failures=False)

    if std_junction.is_standardized:
        return util.safe_return_junction(std_junction.junction, row["aa"]), std_junction.j_gene_match

    # Next: allow non-functional genes (these are being used in the IEDB by some receptors)
    std_junction = tt.junction.standardize(seq=input_cdr3,
                                           locus=get_locus(row['chain_type']),
                                           j_symbol=j_symbol,
                                           v_symbol=v_symbol,
                                           species=species,
                                           allow_c_correction=True if v_symbol is not None else False,
                                           allow_fw_correction=True if j_symbol is not None else False,
                                           max_v_reconstruction=3 if v_symbol is not None and is_curated_cdr3 else 1,
                                           max_j_reconstruction=3 if j_symbol is not None and is_curated_cdr3 else 1,
                                           enforce_functional_v=False,
                                           enforce_functional_j=False,
                                           log_failures=False)

    if std_junction.is_standardized:
        if j_symbol is not None or util.matches_basic_junction_pattern(std_junction.junction):
            return util.safe_return_junction(std_junction.junction, row["aa"]), std_junction.j_gene_match

    if is_curated_cdr3:
        # todo consider dropping this part

        # if the CDR3 is curated (i.e., not retrieved from full aa seq), try with any gene, but allow no error correction/reconstruction
        std_junction_no_genes = tt.junction.standardize(seq=input_cdr3 if std_junction.attempted_fix is None else std_junction.attempted_fix,
                                                        locus=get_locus(row['chain_type']),
                                                        j_symbol=None,
                                                        v_symbol=None,
                                                        species=species,
                                                        allow_c_correction=False,
                                                        allow_fw_correction=False,
                                                        max_v_reconstruction=1,
                                                        max_j_reconstruction=1,
                                                        enforce_functional_v=True,
                                                        enforce_functional_j=True,
                                                        log_failures=False)

        if std_junction_no_genes.is_standardized:
            logging.info(f"Corrected CDR3/Junction without V/J gene info: \"{std_junction_no_genes.junction}\"")
            return util.safe_return_junction(std_junction_no_genes.junction, row["aa"]), std_junction_no_genes.j_gene_match

        if std_junction.attempted_fix is not None:
            if util.matches_basic_junction_pattern(std_junction.attempted_fix) and is_curated_cdr3:
                return util.safe_return_junction(std_junction.attempted_fix, row["aa"]), std_junction.j_gene_match

        logging.warning(f"Failed calculating {row['species_latin']} {row['chain_type']} CDR3/Junction based on Curated CDR3 \"{row["cdr3_seq_curated"]}\", V \"{row["vgene_curated"]}\", J \"{row["jgene_curated"]}\": {std_junction.error}")
        return None, None

def safe_standardize_junction(row, vgene_calc, jgene_calc, species):
    input_cdr3, is_curated_cdr3 = get_input_cdr3(row)

    if input_cdr3 is None:
        return None, None
    try:
        return tt_standardize_junction(row, vgene_calc, jgene_calc, species, input_cdr3, is_curated_cdr3)
    except Exception as e:
        logging.error(
            f"Failed calculating {row['chain_type']} CDR3/Junction based on CDR3 \"{row["cdr3_seq_curated"]}\", V \"{row["vgene_curated"]}\", J \"{row["jgene_curated"]}\": {e}")

    return None, None

def correct_cdr_to_germline(germline_cdrs, curated_cdr, cdr_type, receptor_type):
    if pd.notna(curated_cdr) and len(germline_cdrs) > 0:
        if curated_cdr in germline_cdrs:
            return curated_cdr

        for valid_cdr in sorted(germline_cdrs, key=len)[::-1]:
            if len(valid_cdr) >= 3 and ((valid_cdr in curated_cdr) or (curated_cdr in valid_cdr)):
                logging.info(f"Curated {cdr_type} ({curated_cdr}) has been corrected to valid germline subsequence {valid_cdr}")
                return valid_cdr

        if receptor_type == "TR": # todo: curated 'LIRSNE' germlines {'IRSNERE'} -> if there is only one germline, try overlap?
            logging.warning(f"Curated {cdr_type} {curated_cdr} does not appear in possible germline {cdr_type}s for this gene: {germline_cdrs}. This could mean the V gene or {cdr_type} is incorrect.")


def safe_standardize_cdr12(curated_cdr1, curated_cdr2, module, vgene_calc, species):
    germline_cdr1s, germline_cdr2s = set(), set()

    for v_gene in vgene_calc.split(util.GENE_SEP):
        # Important: do not allow subgroup-level information for standardization (not sufficiently specific)
        std_v_gene = module.standardize(symbol=v_gene, species=species, log_failures=False, allow_subgroup=False, enforce_functional=False)

        if std_v_gene.is_standardized:
            germline_cdr1s.update(std_v_gene.get_aa_sequences("CDR1-IMGT", enforce_functional=True).values())
            germline_cdr2s.update(std_v_gene.get_aa_sequences("CDR2-IMGT", enforce_functional=True).values())

    germline_cdr1s = {cdr for cdr in germline_cdr1s if cdr is not None and cdr != ""}
    germline_cdr2s = {cdr for cdr in germline_cdr2s if cdr is not None and cdr != ""}

    cdr1_standardized = correct_cdr_to_germline(germline_cdr1s, curated_cdr1, "CDR1", vgene_calc[0:2])
    cdr2_standardized = correct_cdr_to_germline(germline_cdr2s, curated_cdr2, "CDR2", vgene_calc[0:2])

    if util.GENE_SEP not in vgene_calc:
        if std_v_gene.is_standardized and std_v_gene.receptor_type == "TR":
            if cdr1_standardized is None and len(germline_cdr1s) == 1:
                cdr1_standardized = germline_cdr1s.pop()
                logging.info(f"No CDR1 available, assigning germline CDR1 {cdr1_standardized} from {std_v_gene}")
            if cdr2_standardized is None and len(germline_cdr2s) == 1:
                cdr2_standardized = germline_cdr2s.pop()
                logging.info(f"No CDR2 available, assigning germline CDR2 {cdr2_standardized} from {std_v_gene}")

    return cdr1_standardized, cdr2_standardized

def standardize_row(row):
    with set_logging_context(template_row=row["row"], chain=row["chain"]):
        module = get_tr_ig_package(row['chain_type'])

        if module is None:
            return None, None, None, None, None, None, None

        species = util.safe_get_field_from_df(row["species_latin"])

        vgene_calc = safe_standardize_gene(module, row["vgene_curated"], species)
        dgene_calc = safe_standardize_gene(module, row["dgene_curated"], species)
        jgene_calc = safe_standardize_gene(module, row["jgene_curated"], species)

        junction_calc, junction_jgene = safe_standardize_junction(row, vgene_calc, jgene_calc, species)

        cdr1_calc, cdr2_calc, cdr3_calc = None, None, None

        if junction_calc:
            if len(junction_calc) > 2:
                cdr3_calc = junction_calc[1:-1]

                if jgene_calc is None and junction_jgene is not None and junction_jgene != "":
                    jgene_calc = junction_jgene
            else:
                logging.warning(f"Calculated junction is too short: '{junction_calc}' (original input: '{row['cdr3_seq_curated']}')")

        # Only validate/retrieve germline CDR1 and CDR2 if they cannot be calculated from full seq with ANARCII
        if pd.isna(row["aa"]) and vgene_calc is not None:
            cdr1_calc, cdr2_calc = safe_standardize_cdr12(row['cdr1_seq_curated'], row['cdr2_seq_curated'], module, vgene_calc, species)

    return vgene_calc, dgene_calc, jgene_calc, cdr1_calc, cdr2_calc, cdr3_calc, junction_calc

def fix_with_tidytcells(full_input):
    full_input = full_input.copy()

    new_columns = ["vgene_calculated", "dgene_calculated", "jgene_calculated", "cdr1_seq_calculated", "cdr2_seq_calculated", "cdr3_seq_calculated", "junction_seq_calculated"]

    full_input[new_columns] = pd.DataFrame(full_input.apply(standardize_row, axis=1).tolist(), index=full_input.index)

    return full_input


