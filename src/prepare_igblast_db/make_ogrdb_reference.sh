if [ -z "$1" ]; then
    echo "Error: Missing argument IGBLAST_FOLDER"
    echo "Usage: $0 path/to/folder"
    exit 1
fi


IGBLAST_ROOT=$1
OGRDB_ROOT=$PWD/igblast_db_ogrdb


mkdir $OGRDB_ROOT
cd $OGRDB_ROOT

mkdir vdj_aa
mkdir database
mkdir optional_file


# Downloading germline sets
# - human
download_germline_set "Homo sapiens" IGH -f MULTI-IGBLAST -n IGH_VDJ
download_germline_set "Homo sapiens" IGL -f MULTI-IGBLAST
download_germline_set "Homo sapiens" IGK -f MULTI-IGBLAST
# - macaque
download_germline_set "Macaca mulatta" IGH -f MULTI-IGBLAST
download_germline_set "Macaca mulatta" IGL -f MULTI-IGBLAST
download_germline_set "Macaca mulatta" IGK -f MULTI-IGBLAST
# - mouse: todo (per strain)
# download_germline_set "Mus musculus" IGH -f MULTI-IGBLAST


# Remove TMP names with *i and _
cat Homo_sapiens_IGH_VDJ_IGH.ndm > Homo_sapiens_IGH_VDJ_IGH_withtmpnames.ndm
cat Homo_sapiens_IGH_VDJ_IGH_V.fasta > Homo_sapiens_IGH_VDJ_IGH_V_withtmpnames.fasta
cat Homo_sapiens_IGH_VDJ_IGH_V_gapped.fasta > Homo_sapiens_IGH_VDJ_IGH_V_gapped_withtmpnames.fasta
cat Homo_sapiens_IGL.ndm >  Homo_sapiens_IGL_withtmpnames.ndm
cat Homo_sapiens_IGL_V.fasta >  Homo_sapiens_IGL_V_withtmpnames.fasta
cat Homo_sapiens_IGL_V_gapped.fasta >  Homo_sapiens_IGL_V_gapped_withtmpnames.fasta

grep -v '\*i' Homo_sapiens_IGH_VDJ_IGH_withtmpnames.ndm > Homo_sapiens_IGH_VDJ_IGH.ndm
awk '/^>/{drop=/\*i/} !drop' Homo_sapiens_IGH_VDJ_IGH_V_withtmpnames.fasta > Homo_sapiens_IGH_VDJ_IGH_V.fasta
awk '/^>/{drop=/\*i/} !drop' Homo_sapiens_IGH_VDJ_IGH_V_gapped_withtmpnames.fasta > Homo_sapiens_IGH_VDJ_IGH_V_gapped.fasta

grep -v '_' Homo_sapiens_IGL_withtmpnames.ndm > Homo_sapiens_IGL.ndm
awk '/^>/{drop=/_/} !drop' Homo_sapiens_IGL_V_withtmpnames.fasta > Homo_sapiens_IGL_V.fasta
awk '/^>/{drop=/_/} !drop' Homo_sapiens_IGL_V_gapped_withtmpnames.fasta > Homo_sapiens_IGL_V_gapped.fasta




download_germline_set "Mus musculus" IGH -f MULTI-IGBLAST -n "BALB/c IGH"
download_germline_set "Mus musculus" IGH -f MULTI-IGBLAST -n "C57BL/6 IGH"
download_germline_set "Mus musculus" IGH -f MULTI-IGBLAST -n "C57BL/6 IGHV"
download_germline_set "Mus musculus" IGH -f MULTI-IGBLAST -n "CAST/EiJ IGH"
download_germline_set "Mus musculus" IGH -f MULTI-IGBLAST -n "LEWES/EiJ IGH"
download_germline_set "Mus musculus" IGH -f MULTI-IGBLAST -n "MSM/MsJ IGH"
download_germline_set "Mus musculus" IGH -f MULTI-IGBLAST -n "NOD/ShiLtJ IGH"
download_germline_set "Mus musculus" IGH -f MULTI-IGBLAST -n "PWD/PhJ IGH"


# BALB/c IGH, C57BL/6 IGH, C57BL/6 IGHV, CAST/EiJ IGH, LEWES/EiJ IGH, MSM/MsJ IGH, NOD/ShiLtJ IGH, PWD/PhJ IGH
# 129S1/SvlmJ, A/J IGLV, AKR/J IGLV, BALB/cByJ IGLV, C3H/HeJ IGLV, C57BL/6J IGLV, CAST/EiJ IGLV, CBA/J IGLV, DBA/1J IGLV, DBA/2J IGLV, LEWES/EiJ IGLV, MRL/MpJ IGLV, MSM/MsJ IGLV, NOD/ShiLtJ IGLV, NOR/LtJ IGLV, NZB/BlNJ IGLV, PWD/PhJ IGLV, SJL/J IGLV, 129S1/SvlmJ IGLV, IGLJ (all strains), BALB/c/ByJ IGLV, 129S1/SvImJ IGLV
# 129S1/SvImJ IGKV, AKR/J IGKV, BALB/c/ByJ IGKV, C3H/HeJ IGKV, C57BL/6J IGKV, CAST/EiJ IGKV, CBA/J IGKV, DBA/1J IGKV, DBA/2J IGKV, LEWES/EiJ IGKV, MRL/MpJ IGKV, MSM/MsJ IGKV, NOD/ShiLtJ IGKV, NOR/LtJ IGKV, NZB/BlNJ IGKV, PWD/PhJ IGKV, SJL/J IGKV, A/J IGKV, IGKJ (all strains)

# Combine all IG V, IG J (to discover best matching locus)
# - human
cat Homo_sapiens_IGH_VDJ_IGH_V.fasta Homo_sapiens_IGL_V.fasta Homo_sapiens_IGK_V.fasta > combined_Homo_sapiens_IG_V.fasta
cat Homo_sapiens_IGH_VDJ_IGH_J.fasta Homo_sapiens_IGL_J.fasta Homo_sapiens_IGK_J.fasta > combined_Homo_sapiens_IG_J.fasta
# - macaque
cat Macaca_mulatta_IGH_V.fasta Macaca_mulatta_IGL_V.fasta Macaca_mulatta_IGK_V.fasta > combined_Macaca_mulatta_IG_V.fasta
cat Macaca_mulatta_IGH_J.fasta Macaca_mulatta_IGL_J.fasta Macaca_mulatta_IGK_J.fasta > combined_Macaca_mulatta_IG_J.fasta



# Combining ndm files for translation
cat Homo_sapiens_IGH_VDJ_IGH.ndm <(tail -n +2 Homo_sapiens_IGK.ndm) <(tail -n +2 Homo_sapiens_IGL.ndm) > combined_Homo_sapiens_IG.ndm
cat Macaca_mulatta_IGH.ndm <(tail -n +2 Macaca_mulatta_IGK.ndm) <(tail -n +2 Macaca_mulatta_IGL.ndm) > combined_Macaca_mulatta_IG.ndm



# translating the ungapped fasta files to make blast database
# - human
python $ROOT/src/prepare_igblast_db/translate_ogrdb.py -i $OGRDB_ROOT/combined_Homo_sapiens_IG_V.fasta -o $OGRDB_ROOT/combined_aa_Homo_sapiens_IG_V.fasta
# - macaque
python $ROOT/src/prepare_igblast_db/translate_ogrdb.py -i $OGRDB_ROOT/combined_Macaca_mulatta_IG_V.fasta -o $OGRDB_ROOT/combined_aa_Macaca_mulatta_IG_V.fasta



# Making 'optional files' -> J gene start indices
cat Homo_sapiens_IGH_VDJ_IGH.aux <(tail -n +2 Homo_sapiens_IGK.aux) <(tail -n +2 Homo_sapiens_IGL.aux) > optional_file/human_gl.aux
cat Macaca_mulatta_IGH.aux <(tail -n +2 Macaca_mulatta_IGK.aux) <(tail -n +2 Macaca_mulatta_IGL.aux) > optional_file/rhesus_monkey_gl.aux



# Making IgBLAST nt databases
# - human
$IGBLAST_ROOT/bin/makeblastdb -parse_seqids -dbtype nucl -in $OGRDB_ROOT/combined_Homo_sapiens_IG_V.fasta -out $OGRDB_ROOT/database/ogrdb_human_ig_v
$IGBLAST_ROOT/bin/makeblastdb -parse_seqids -dbtype nucl -in $OGRDB_ROOT/Homo_sapiens_IGH_VDJ_IGH_D.fasta -out $OGRDB_ROOT/database/ogrdb_human_ig_d
$IGBLAST_ROOT/bin/makeblastdb -parse_seqids -dbtype nucl -in $OGRDB_ROOT/combined_Homo_sapiens_IG_J.fasta -out $OGRDB_ROOT/database/ogrdb_human_ig_j
# - macaque
$IGBLAST_ROOT/bin/makeblastdb -parse_seqids -dbtype nucl -in $OGRDB_ROOT/combined_Macaca_mulatta_IG_V.fasta -out $OGRDB_ROOT/database/ogrdb_rhesus_monkey_ig_v
$IGBLAST_ROOT/bin/makeblastdb -parse_seqids -dbtype nucl -in $OGRDB_ROOT/Macaca_mulatta_IGH_D.fasta -out $OGRDB_ROOT/database/ogrdb_rhesus_monkey_ig_d
$IGBLAST_ROOT/bin/makeblastdb -parse_seqids -dbtype nucl -in $OGRDB_ROOT/combined_Macaca_mulatta_IG_J.fasta -out $OGRDB_ROOT/database/ogrdb_rhesus_monkey_ig_j


# Making IgBLAST aa databases
# - human
$IGBLAST_ROOT/bin/makeblastdb -parse_seqids -dbtype prot -in $OGRDB_ROOT/combined_aa_Homo_sapiens_IG_V.fasta -out $OGRDB_ROOT/database/ogrdb_aa_human_ig_v
# - macaque
$IGBLAST_ROOT/bin/makeblastdb -parse_seqids -dbtype prot -in $OGRDB_ROOT/combined_aa_Macaca_mulatta_IG_V.fasta -out $OGRDB_ROOT/database/ogrdb_aa_rhesus_monkey_ig_v



# translating the gapped fasta files for MakeDB.py input
# - human
python $ROOT/src/prepare_igblast_db/translate_ogrdb.py -i $OGRDB_ROOT/Homo_sapiens_IGH_VDJ_IGH_V_gapped.fasta -o $OGRDB_ROOT/vdj_aa/ogrdb_aa_human_IGHV.fasta
python $ROOT/src/prepare_igblast_db/translate_ogrdb.py -i $OGRDB_ROOT/Homo_sapiens_IGK_V_gapped.fasta -o $OGRDB_ROOT/vdj_aa/ogrdb_aa_human_IGKV.fasta
python $ROOT/src/prepare_igblast_db/translate_ogrdb.py -i $OGRDB_ROOT/Homo_sapiens_IGL_V_gapped.fasta -o $OGRDB_ROOT/vdj_aa/ogrdb_aa_human_IGLV.fasta
# - macaque
python $ROOT/src/prepare_igblast_db/translate_ogrdb.py -i $OGRDB_ROOT/Macaca_mulatta_IGH_V_gapped.fasta -o $OGRDB_ROOT/vdj_aa/ogrdb_aa_rhesus_monkey_IGHV.fasta
python $ROOT/src/prepare_igblast_db/translate_ogrdb.py -i $OGRDB_ROOT/Macaca_mulatta_IGK_V_gapped.fasta -o $OGRDB_ROOT/vdj_aa/ogrdb_aa_rhesus_monkey_IGKV.fasta
python $ROOT/src/prepare_igblast_db/translate_ogrdb.py -i $OGRDB_ROOT/Macaca_mulatta_IGL_V_gapped.fasta -o $OGRDB_ROOT/vdj_aa/ogrdb_aa_rhesus_monkey_IGLV.fasta

