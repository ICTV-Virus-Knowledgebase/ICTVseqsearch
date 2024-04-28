#!/usr/bin/env bash
#
# extract GenBank accessions from VMR.xlsx into a text file
# query NCBI for all the matching RefSeq accessions
#
VMR="VMR_MSL39_v1.xlsx"
VMR_E="processed_accessions_e.tsv"
ACC_E="processed_accessions_e_acc_only.txt"
OUT_DIR="refseq"
PREFIX="$OUT_DIR/e_genbank.txt."
PREFIX_OUT="$OUT_DIR/e_refseq.txt."
OUT_MAP="genbank_refseq_map.txt"
TEST1="test_acc_list_1.txt"
TEST10="test_acc_list_10.txt"

if [ "$1" == "-h" ]; then
    echo "SYNTAX: $0 [-h] [-query]"
    echo " "
    echo " -query: only query NCBI for missing "
    exit 0
fi

QUERY_ONLY=0
if [ "$1" == "-query" ]; then
    QUERY_ONLY=1
fi

##
## PARSE MSL, extract genbank accession list
##
if [ ! $QUERY_ONLY ]; then
    
    conda activate  conda/vmr_openpyxl3

    if [[ $VMR -nt $VMR_E ]]; then 
	echo RUN ./VMR_to_fasta.py -verbose -VMR_file_name $VMR -mode VMR -ea E
	./VMR_to_fasta.py -verbose -VMR_file_name $VMR -mode VMR -ea E
    else
	echo "SKIP: $VMR_E up-to-date"
    fi
    wc -l $VMR_E

    # just accessions
    cut -f 3   $VMR_E > $ACC_E
    wc -l $ACC_E

    # split main file into small piecees
    mkdir -p $OUT_DIR
    split -l 100  $ACC_E $PREFIX

    # short test files
    head -2  $ACC_E  | grep -v Accession_IDs > $TEST1
    head -11  $ACC_E | grep -v Accession_IDs > $TEST10

fi

##
## QUERY NCBI
##

 
# NCBI query for each file
mkdir -p test
for file in $TEST1 $TEST10 ; do
    out=test/${file}.out
    echo "############### $file -> $out ##############"
    if [[ ${file} -nt ${out} ]]; then 
	echo "QUERY NCBI:"
	elink -input $file -db nucleotide -target nucleotide -name nuccore_nuccore_gbrs \
	    | efetch -format docsum \
	    | xtract -pattern DocumentSummary -sep "\t" -element Caption AssemblyAcc Organism TaxID \
	    | tee $out
    else
	echo "SKIP"
    fi
done



# NCBI query for each file
for file in ${PREFIX}* ; do
    out=$(echo $file | perl -pe "s|$PREFIX|$PREFIX_OUT|;")
    echo -n "# [$file => $out] ##  "
    if [[ ! -e ${out} ]]; then
	echo -n "QUERY: "
	elink -input $file -db nucleotide -target nucleotide -name nuccore_nuccore_gbrs \
	    | efetch -format docsum \
	    | xtract -pattern DocumentSummary -sep "\t" -element Caption AssemblyAcc Organism TaxID \
		     2>&1 > $out
	wc -l $out
    else
	echo "SKIP"
    fi
done

#
# merge output files
#
sort -k2,2 ${PREFIX_OUT}* > $OUT_MAP
echo "FROM $VMR"
echo "   GENBANK: $(wc -l $VMR_E)"
echo "   REFSEQ:  $(wc -l $OUT_MAP)"
