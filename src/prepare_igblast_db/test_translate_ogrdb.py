from prepare_igblast_db.translate_ogrdb import *


def test_trim_to_frame():
    assert trim_to_frame("AAACCCTTT", frame_start=0) == "AAACCCTTT"
    assert trim_to_frame("AAACCCTTT", frame_start=1) == "AACCCT"
    assert trim_to_frame("AAACCCTTT", frame_start=2) == "ACCCTT"


def test_translate_genes():
    genes_nt = SequenceEntry(name=["TRAV1-1", "TRAV2-1", "TRAV3-1"],
                             sequence=["AAACCCTTTG","ACCCTTTG","ACCCTTTGTGA"],)


    frame_start_dict = {"TRAV1-1":0, "TRAV2-1":1, "TRAV3-1":2}

    output = translate_genes(genes_nt, frame_start_dict)
    
    assert set([str(name) for name in output.name]) == set(frame_start_dict.keys())
    assert str(output.sequence[0]) == "KPF"
    assert str(output.sequence[1]) == "PF"
    assert str(output.sequence[2]) == "PL*"
