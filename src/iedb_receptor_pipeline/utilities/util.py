import os
from pathlib import Path

import logging
import pandas as pd
import shutil
import re


PROJECT_ROOT = None
IMGT_GERMLINE_ROOT = None
OGRDB_GERMLINE_ROOT = None


IGBLAST_ROOT = PROJECT_ROOT / "IgBLAST/ncbi-igblast-1.22.0"
IGBLAST_N = IGBLAST_ROOT / "bin/igblastn"
IGBLAST_P = IGBLAST_ROOT / "bin/igblastp"

NT_ALPHABET = {'g', 't', 'a', 'c'}
NT_ALPHABET_AMBIGUOUS = NT_ALPHABET.union({'n'})
AA_ALPHABET_STRICT = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
AA_ALPHABET_AMBIGUOUS = AA_ALPHABET_STRICT.union({'X', '*'})


# see also: https://dst.lji.org/Organism.html

NCBITAXON_TO_LATIN = {"9606": "Homo sapiens",
                     "10090": "Mus musculus",
                     "30538": "Vicugna pacos",
                     "9986": "Oryctolagus cuniculus",
                     "9844": "Lama glama",
                     "9541": "Macaca fascicularis",
                     "9838": "Camelus dromedarius",
                     "10116": "Rattus norvegicus",
                     "9544": "Macaca mulatta",
                     "10114": "Rattus",
                     "9837": "Camelus bactrianus",
                      "36176": "Chiloscyllium plagiosum",
                      "9835": "Camelidae",
                      "9031": "Gallus gallus",
                      "9598": "Pan troglodytes",
                      "9839": "Lama",
                      "7801": "Ginglymostoma cirratum",
                      "9913": "Bos taurus",
                      "9823": "Sus scrofa",
                      "9534": "Chlorocebus aethiops",
                      "7797": "Squalus acanthias",
                      "9980": "Lepus",
                      "168098": "Orectolobus maculatus",
                      "10088": "Mus",
                      "447135": "Myodes glareolus",
                      "419612": "Camelus ferus",
                      "3122392": "Nothocricetulus migratorius",
                      "32630": "synthetic construct"
                      }

LATIN_TO_NCBITAXON = {v: k for k, v in NCBITAXON_TO_LATIN.items()}


ORGANISM_ID_TO_LATIN = {
    "10002763": "Mus musculus",  # Mus musculus HLA-A*11:01 Tg
    "10000616": "Mus musculus",  # Mus musculus HLA-DRB1*04:01 Tg
    "10003310": "Mus musculus",  # Mus musculus HLA-A*03:01 Tg
    "477816": "Mus musculus",
    "10001540": "Mus musculus",
    "1385377": "Mus musculus",
    "10000649": "Mus musculus",
    "10000645": "Mus musculus",
    "10000973": "Mus musculus",
    "10000028": "Mus musculus",
    "10001241": "Mus musculus",
    "10000053": "Mus musculus",
    "10000632": "Mus musculus",
    "10002278": "Mus musculus",
    "10000067": "Mus musculus",
    "10000653": "Mus musculus",
    "10001537": "Mus musculus",
    "10001538": "Mus musculus",
    "10000061": "Mus musculus",
    "10000000": "Mus musculus",
    "10000253": "Mus musculus",
    "10000228": "Mus musculus",
    "10002138": "Mus musculus",
    "10002377": "Mus musculus",
    "10002378": "Mus musculus",
    "10002379": "Mus musculus",
    "10001539": "Mus musculus",
    "10091": "Mus musculus",
    "10092": "Mus musculus",
    "39442": "Mus musculus",
    "35531": "Mus musculus",
    "46456": "Mus musculus",
    "57486": "Mus musculus",
    "80274": "Mus musculus",
    "116058": "Mus musculus",
    "179238": "Mus musculus",
    "1266728": "Mus musculus",
    "3109569": "Mus musculus",
    "477815": "Mus musculus",
    "947985": "Mus musculus",
    "1643390": "Mus musculus",
    "1879032": "Mus musculus",
    "10000054": "Mus musculus",
    "10000063": "Mus musculus",
    "10000106": "Mus musculus",
    "10000206": "Oryctolagus cuniculus",
    "10000230": "Mus musculus",
    "10000233": "Mus musculus",
    "10000236": "Mus musculus",
    "10000263": "Mus musculus",
    "10000274": "Mus musculus",
    "10000289": "Mus musculus",
    "10000635": "Mus musculus",
    "10001019": "Mus musculus",
    "10001758": "Mus musculus"
}



ONTIE_TO_LATIN = {
    "0000001": "Mus musculus",
    "0000037": "Mus musculus",
    "0000050": "Mus musculus",
    "0000056": "Mus musculus",
    "0000058": "Mus musculus",
    "0000062": "Mus musculus",
    "0000147": "Mus musculus",
    "0000149": "Mus musculus",
    "0000152": "Mus musculus",
    "0000155": "Mus musculus",
    "0000163": "Mus musculus",
    "0000178": "Mus musculus",
    "0000189": "Mus musculus",
    "0000190": "Mus musculus",
    "0000540": "Mus musculus",
    "0000543": "Mus musculus",
    "0000552": "Mus musculus",
    "0000556": "Mus musculus",
    "0000561": "Mus musculus",
    "0000884": "Mus musculus",
    "0001343": "Mus musculus",
    "0000126": "Oryctolagus cuniculus",
    "0000573": "Rattus norvegicus",
    "0000525": "Mus musculus", # HuMAb Medarex",
    "0000563": "Mus musculus", # Xenomouse",
    "0000807": "Mus musculus", # Kunming",
    "0001345": "Mus musculus", # TCR transgenic",
    "0001558": "Mus musculus", # 6218 TCR Tg",
    "0000101": "Mus musculus", # HY-A1 TCR beta chain Tg",
    "0001408": "Mus musculus", # RBF",
    "0000204": "Mus musculus", # BALB.B",
}

CURATION_TEMPLATE_CURATED_COLUMNS = ["curated_receptor_id", "ASSAY_ID", "REFERENCE_ID", "ref_name", "synonyms", "RR_ID", "receptor_type",
                                     "curated_chain1_id", "chain1_species", "chain1_type", "chain1_nucleotide", "chain1_nt_accession", "chain1_v_gene_curated", "chain1_d_gene_curated", "chain1_j_gene_curated", "chain1_pro_seq", "chain1_accession", "chain1_cdr1_seq_curated", "chain1_cdr2_seq_curated", "chain1_cdr3_seq_curated", "chain1_vdomain_curated",
                                     "curated_chain2_id", "chain2_species", "chain2_type", "chain2_nucleotide", "chain2_nt_accession", "chain2_v_gene_curated", "chain2_d_gene_curated", "chain2_j_gene_curated", "chain2_pro_seq", "chain2_accession", "chain2_cdr1_seq_curated", "chain2_cdr2_seq_curated", "chain2_cdr3_seq_curated", "chain2_vdomain_curated"]

CURATION_TEMPLATE_CALCULATED_COLUMNS = ["calculated_chain1_id", "chain1_chain_type_calculated", "chain1_v_gene_calculated", "chain1_d_gene_calculated", "chain1_j_gene_calculated", "chain1_cdr1_seq_calculated", "chain1_cdr2_seq_calculated", "chain1_cdr3_seq_calculated", "chain1_junction_seq_calculated", "chain1_vdomain_calculated",
                                        "calculated_chain2_id", "chain2_chain_type_calculated", "chain2_v_gene_calculated", "chain2_d_gene_calculated", "chain2_j_gene_calculated", "chain2_cdr1_seq_calculated", "chain2_cdr2_seq_calculated", "chain2_cdr3_seq_calculated", "chain2_junction_seq_calculated", "chain2_vdomain_calculated"]

CURATION_TEMPLATE_COLUMNS = CURATION_TEMPLATE_CURATED_COLUMNS + CURATION_TEMPLATE_CALCULATED_COLUMNS


LATIN_TO_ANARCI = {"Homo sapiens": "human", "Mus musculus": "mouse", "Mus musculus domesticus": "mouse", "Vicugna pacos": "alpaca", "Oryctolagus cuniculus": "rabbit",
                   "Rattus rattus": "rat", "Macaca mulatta": "rhesus", "Sus domesticus": "pig", "Bos taurus": "cow"}

LATIN_TO_COMMON = {"Homo sapiens": "human", "Mus musculus": "mouse", "Macaca mulatta": "rhesus_monkey"}
COMMON_TO_UNDERSC_LATIN = {"human": "Homo_sapiens",  "mouse": "Mus_musculus"}

BCR_CHAINS = {"heavy", "light", "kappa_light", "lambda_light"}
TCR_CHAINS = {"alpha", "beta", "gamma", "delta"}

RECEPTOR_TYPE_TO_VALID_CHAINS = {"heavylight": BCR_CHAINS,
                                 "heavyheavy": {"heavy"},
                                 "heavy": {"heavy"},
                                 "light": {"light", "kappa_light", "lambda_light"},
                                 "lightlight": {"light", "kappa_light", "lambda_light"},
                                 "alphabeta": {"alpha", "beta"},
                                 "gammadelta": {"gamma", "delta"},
                                 "scfv": BCR_CHAINS,
                                 "tscfv": TCR_CHAINS,
                                 "construct": {*{"construct"}, *BCR_CHAINS, *TCR_CHAINS},
                                 "ignar": {"ignar", "heavy"}}


TOOL_INPUT_COLUMNS = ['curated_receptor_id', 'row', 'chain_id_cur', 'chain_id_calc', 'receptor_type', 'chain',
                        'species_latin', 'organism_name', 'aa', 'nt', 'chain_type', 'v_dom_seq', 'vgene_curated',
                        'dgene_curated', 'jgene_curated', 'cdr1_seq_curated', 'cdr2_seq_curated', 'cdr3_seq_curated']


CHAIN_LOCUS_TO_IEDB = {"TRA": "alpha", "TRB": "beta", "TRG": "gamma", "TRD": "delta", "IGH": "heavy", "IGL": "lambda_light", "IGK": "kappa_light"}

CHAIN_IEDB_TO_ANARCI = {"alpha": ["A"], "beta": ["B"], "gamma": ["G"], "delta": ["D"], "heavy": ["H"], "light": ["L", "K"], "kappa_light": ["K"], "lambda_light": ["L"]}
CHAIN_ANARCI_TO_IEDB_ALL = {"H": ["heavy"], "L": ["lambda_light", "light"], "A": ["alpha"], "B": ["beta"], "K": ["kappa_light", "light"], "D": ["delta"], "G": ["gamma"]}
CHAIN_ANARCI_TO_IEDB_SINGLE = {"H": "heavy", "L": "lambda_light", "A": "alpha", "B": "beta", "K": "kappa_light", "D": "delta", "G": "gamma"}

CHAIN_TO_LEGAL_LOCI = {"alpha": {"TRA", "TRD"}, # , "TRD"
                       "beta": {"TRB"},
                       "gamma": {"TRG"},
                       "delta": {"TRD"}, # , "TRA"
                       "heavy": {"IGH"},
                       "light": {"IGL", "IGK"},
                       "kappa_light": {"IGK"},
                       "lambda_light": {"IGL"}}


VDJ_CHAINS = {"beta", "delta", "heavy"}
VJ_CHAINS = {"alpha", "gamma", "light", "kappa_light", "lambda_light"}


IEDB_SHARK_SPECIES = ["Ginglymostoma cirratum", "Orectolobus macalatus", "Squalus acanthias", "Chiloscyllium plagiosum"]
IEDB_CAMEL_SPECIES = ["Lama glama", "Camelus dromedarius", "Camelidae", "Vicugna pacos", "Lama", "Camelus bactrianus", "Camelus ferus"]


EVAL_TRUST_CUTOFF = 1e-06

GENE_SEP = ", "

# CDR3 range(105, 118) -> means IMGT positions 105 until 117 are included (118 excluded)
CDR_RANGES = {"cdr1": range(27, 39), "cdr2": range(56, 66), "cdr3": range(105, 118), "junction": range(104, 119)}


DOUBLE_TO_SINGLE_HEADER = {('Receptor', 'Group IRI'): 'receptor_group_id',
                           ('Receptor', 'IEDB Receptor ID'): 'receptor_iedb_id',
                           ('Receptor', 'Reference Name'): 'REFERENCE_ID',
                           ('Receptor', 'Type'): 'receptor_type',
                           ('Reference', 'IEDB IRI'): 'reference_iri',
                           ('Epitope', 'IEDB IRI'): 'epitope_iedb_id',
                           ('Epitope', 'Name'): 'epitope_name',
                           ('Epitope', 'Source Molecule'): 'epitope_source_molecule',
                           ('Epitope', 'Source Organism'): 'epitope_source_organism',
                           ('Assay', 'Type'): 'assay_type',
                           ('Assay', 'IEDB IDs'): 'assay_iedb_ids',
                           ('Assay', 'MHC Allele Names'): 'mhc_allele_names',
                           ('Chain 1', 'Type'): 'chain1_type',
                           ('Chain 1', 'Organism IRI'): 'chain1_species',
                           ('Chain 1', 'Nucleotide Sequence'): 'chain1_nucleotide',
                           ('Chain 1', 'Curated V Gene'): 'chain1_v_gene_curated',
                           ('Chain 1', 'Curated D Gene'): 'chain1_d_gene_curated',
                           ('Chain 1', 'Curated J Gene'): 'chain1_j_gene_curated',
                           ('Chain 1', 'Protein Sequence'): 'chain1_pro_seq',
                           ('Chain 1', 'Protein IRI'): 'chain1_accession',
                           ('Chain 1', 'CDR3 Curated'):  'chain1_cdr3_seq_curated',
                           ('Chain 1', 'CDR1 Curated'): 'chain1_cdr1_seq_curated',
                           ('Chain 1', 'CDR2 Curated'): 'chain1_cdr2_seq_curated',
                           ('Chain 2', 'Type'): 'chain2_type',
                           ('Chain 2', 'Organism IRI'): 'chain2_species',
                           ('Chain 2', 'Nucleotide Sequence'): 'chain2_nucleotide',
                           ('Chain 2', 'Curated V Gene'): 'chain2_v_gene_curated',
                           ('Chain 2', 'Curated D Gene'): 'chain2_d_gene_curated',
                           ('Chain 2', 'Curated J Gene'): 'chain2_j_gene_curated',
                           ('Chain 2', 'Protein Sequence'): 'chain2_pro_seq',
                           ('Chain 2', 'Protein IRI'): 'chain2_accession',
                           ('Chain 2', 'CDR3 Curated'): 'chain2_cdr3_seq_curated',
                           ('Chain 2', 'CDR1 Curated'): 'chain2_cdr1_seq_curated',
                           ('Chain 2', 'CDR2 Curated'): 'chain2_cdr2_seq_curated'}


def build_path(path, delete_if_exists=False):
    path = Path(path)

    if delete_if_exists:
        delete_file_or_folder(path)

    if "." in path.name:
        path = path.parent

    if not path.is_dir():
        os.makedirs(path)

def delete_file_or_folder(path):
    path = Path(path)

    if path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)

def str_format_unique_list(input_iterable):
    return ", ".join([str(item) for item in sorted({str(item) for item in input_iterable})])

def get_alleles_from_fasta_header(file_location):
    with open(file_location, "r") as f:
        headers = {line[1:].strip() for line in f.readlines() if line.startswith(">")}
    return headers

def safe_read_csv(path):
    if path.exists():
        return pd.read_csv(path)

def safe_get_field_from_df(field):
    if pd.notna(field):
        return field

def remove_allele(field):
    return field.rsplit("*")[0]

def split_scfv(seq, linker_len=8):
    if pd.isna(seq):
        return [None, None]

    sg_linker_match = max(re.finditer(r'[SG]{' + str(linker_len) + ',}', seq), key=lambda m: len(m.group()), default=None)

    if sg_linker_match:
        linker = sg_linker_match.group()

        chain1, chain2 = seq.split(linker, maxsplit=1)

        if len(chain1) > 70 and len(chain2) > 70:
            return [chain1, chain2]

    return [None, None]

def is_scfv(seq, linker_len=8):
    chain1, chain2 = split_scfv(seq, linker_len)

    if type(chain1) is str and type(chain2) is str:
        return True

    return False


def is_valid_alphabet(seq, alphabet):
    return all([char in alphabet for char in seq])


def is_valid_junction_in_aa(junction, aa):
    if junction in aa:
        return True

    if junction[-1] in ("W", "F") and aa.endswith(junction[:-1]):
        return True

    logging.info(f"Junction {junction} not in aa {aa}")
    return False

def safe_return_junction(junction_seq, aa_seq):
    if pd.isna(aa_seq):
        return junction_seq

    if is_valid_junction_in_aa(junction_seq, aa_seq):
        return junction_seq

    return None


def matches_basic_junction_pattern(junction):
    if re.match(r"^C.{4,40}[WF]$", junction):
        return True
    return False
