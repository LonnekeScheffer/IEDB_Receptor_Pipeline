from argparse import ArgumentParser

import Bio.Data.CodonTable

from pathlib import Path
import bionumpy as bnp
import pandas as pd


from bionumpy import SequenceEntry


CODON_TABLE_WITH_GAP = {**{"...": ".", "TAA": "*", "TAG": "*", "TGA": "*"},
                        **Bio.Data.CodonTable.standard_dna_table.forward_table}

def parse_arguments():
    args = ArgumentParser()
    args.add_argument("-i", f'--input_file', required=True)
    args.add_argument("-o", f'--output_file', required=True)

    return args.parse_args()

def get_frame_start_per_gene(file, frame_col_name="frame start"):
    assert file.is_file(), f"Expected file {file} does not exist"

    df = pd.read_csv(file, sep="\t", usecols=["#name", frame_col_name])
    return df.set_index("#name")[frame_col_name].to_dict()

def trim_to_frame(sequence_nt, frame_start=0):
    trim_end = len(sequence_nt[frame_start:]) % 3

    if trim_end > 0:
        trimmed = sequence_nt[frame_start:-trim_end]
    else:
        trimmed = sequence_nt[frame_start:]

    assert len(trimmed) > 0, f"Something went wrong while trimming {sequence_nt}"
    return trimmed

def translate_nt_to_aa_with_gaps(dna_seq):
    assert len(dna_seq) % 3 == 0
    return "".join([CODON_TABLE_WITH_GAP[dna_seq[i:i+3]] for i in range(0, len(dna_seq), 3)])

def translate_genes(genes_nt, frame_start_dict):
    trimmed = [trim_to_frame(gene.sequence, frame_start=frame_start_dict[str(gene.name)]) for gene in genes_nt]
    translated = [translate_nt_to_aa_with_gaps(str(seq)) for seq in trimmed]

    return SequenceEntry(name=genes_nt.name, sequence=bnp.as_encoded_array(translated))


def main(args):
    input_file = Path(args.input_file)
    output_file = Path(args.output_file)

    assert input_file.suffix == ".fasta", f"Expected Fasta input file: {input_file}"
    assert output_file.suffix == ".fasta", f"Expected Fasta output file: {output_file}"

    if "Homo_sapiens" in input_file.stem:
        assert "human" in output_file.stem or "Homo_sapiens" in output_file.stem, f"Expected 'human' or 'Homo sapiens' in output file name: {output_file}"
    elif "Macaca_mulatta" in input_file.stem:
        assert "rhesus_monkey" in output_file.stem or "Macaca_mulatta" in output_file.stem, f"Expected 'rhesus_monkey' or 'Macaca_mulatta' in output file name: {output_file}"
    else:
        assert False, f"Supporting 'Homo_sapiens' and 'Macaca_mulatta' only, found neither in input filename: {input_file}"

    if "IGH" in input_file.stem:
        assert 'IGH' in output_file.stem, f"Expected 'IGH' in output file name: {output_file}"
    elif "IGK" in input_file.stem:
        assert "IGK" in output_file.stem, f"Expected 'IGK' in output file name: {output_file}"
    elif "IGL" in input_file.stem:
        assert "IGL" in output_file.stem, f"Expected 'IGL' in output file name: {output_file}"
    elif "IG" in input_file.stem:
        assert "IG" in output_file.stem, f"Expected 'IG' in output file name: {output_file}"
    else:
        assert False, f"Supporting 'IGH', 'IGK', 'IGL', 'IG' only, found neither in input filename: {input_file}"

    assert input_file.is_file(), f"Not a file: {input_file}"
    assert not output_file.is_file(), f"File already exists: {output_file}"

    if input_file.stem.endswith("_gapped"):
        frame_start_file_stem, gene, _ = input_file.stem.rsplit("_", maxsplit=2)
    else:
        frame_start_file_stem, gene = input_file.stem.rsplit("_", maxsplit=1)

    assert gene in ("V", "J"), f"Expected gene to be V or J: {args.gene}"

    if gene == "V":
        frame_start_file = input_file.parent / f"{frame_start_file_stem}.ndm"
        assert frame_start_file.is_file(), f"Not a file: {frame_start_file}"

        frame_start_dict = get_frame_start_per_gene(frame_start_file, frame_col_name="frame start")
    elif gene == "J":
        frame_start_file = input_file.parent / f"{frame_start_file_stem}.aux"
        assert frame_start_file.is_file(), f"Not a file: {frame_start_file}"

        frame_start_dict = get_frame_start_per_gene(frame_start_file, frame_col_name="j_codon_frame")

    genes_nt = bnp.open(input_file).read()

    translated = translate_genes(genes_nt, frame_start_dict)

    with bnp.open(output_file, "w") as f:
        f.write(translated)


if __name__ == '__main__':
    args = parse_arguments()
    main(args)