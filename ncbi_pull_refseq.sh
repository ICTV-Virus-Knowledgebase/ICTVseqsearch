#!/usr/bin/env bash
#
# extract GenBank accessions from VMR.xlsx into a text file
# query NCBI for all the matching RefSeq accessions
#
VMR="VMRs/VMR_MSL39_v1.xlsx"
VMR_TSV="processed_accessions_b.tsv"
ACC_TSV="processed_accessions_b_acc_only.txt"
OUT_DIR="refseq"
PREFIX="$OUT_DIR/b_genbank.txt."
PREFIX_OUT="$OUT_DIR/b_refseq.txt."
OUT_MAP="genbank_refseq_map.txt"
TEST0="test_OQ721911.txt"
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
if [ $QUERY_ONLY -lt 1 ]; then
    
    echo conda activate  conda/vmr_openpyxl3
    conda activate  conda/vmr_openpyxl3

    if [[ $VMR -nt $VMR_TSV ]]; then 
	echo RUN ./VMR_to_fasta.py -verbose -VMR_file_name $VMR -mode VMR -ea B
	./VMR_to_fasta.py -verbose -VMR_file_name $VMR -mode VMR -ea B
    else
	echo "SKIP: $VMR_TSV up-to-date"
    fi
    wc -l $VMR_TSV

    # just accessions - remove ()'s
    cut -f 4   $VMR_TSV \
	| perl -pe 's/\([0-9.]+\)//g;' \
	> $ACC_TSV
    
    wc -l $ACC_TSV

    # split main file into small piecees
    mkdir -p $OUT_DIR
    split -l 100  $ACC_TSV $PREFIX

    # short test files
    head -2  $ACC_TSV  | grep -v Accession_IDs > $TEST1
    head -11  $ACC_TSV | grep -v Accession_IDs > $TEST10

fi

##
## QUERY NCBI
##

 
# NCBI query for each file
mkdir -p test
for file in $TEST0 $TEST1 $TEST10 ; do
    out=test/${file}.out
    echo "############### $file -> $out ##############"
    if [[ ${file} -nt ${out} ]]; then 
	echo "QUERY NCBI:"
	t1=test/${file}.tmp.1.txt
	t2=test/${file}.tmp.2.txt
	t3=$out
	( \
	  elink -input $file -db nucleotide -target nucleotide -name nuccore_nuccore_gbrs > $t1 \
              || echo "ERROR $? : elink -input $file -db nucleotide -target nucleotide -name nuccore_nuccore_gbrs > $t1" > /dev/stderr \
	) \
	    && ( \
		 efetch -format docsum < $t1 > $t2 \
		     || echo "ERROR $? : effect -format docsum < $t1 > $t2"  > /dev/stderr
	) \
	    && ( \
		 cat $t2 \
		     | xtract -pattern DocumentSummary -sep "\t" -element Caption AssemblyAcc Organism TaxID > $t3 \
		     || echo "ERROR $? : xtract -pattern DocumentSummary -sep '\t' -element Caption AssemblyAcc Organism TaxID < $t2 > $t3">/dev/stderr \
	) \
	    && cat $t3
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
echo "   GENBANK: $(wc -l $VMR_TSV)"
echo "   REFSEQ:  $(wc -l $OUT_MAP)"


#merge RefSeq into VMR export
./update_vrm_refseq.py                                                   

# QC check
 cut -f 3 processed_accessions_e.tsv | sort | uniq -c | sort -k1,1n | tail
