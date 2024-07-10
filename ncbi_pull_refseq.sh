#!/usr/bin/env bash
#
# 1. EXTRACT GenBank accessions from VMR.xlsx into a text file
# 2. QUERY: NCBI for all the matching RefSeq accessions
#
#SBATCH --job-name=VMR_NCBI_efetch_refseq_accessions
#SBATCH --output=logs/log.%J.%x.out
#SBATCH --error=logs/log.%J.%x.out
#
# Number of tasks needed for this job. Generally, used with MPI jobs
# Time format = HH:MM:SS, DD-HH:MM:SS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=amd-hdr100 --time=00-12:00:00
##SBATCH --partition=amd-hdr100 --time=06-06:00:00
##SBATCH --partition=medium --time=40:00:00
#
# Number of CPUs allocated to each task. 
#
# Mimimum memory required per allocated  CPU  in  MegaBytes. 
#  last run was 402M
#SBATCH --mem-per-cpu=1000
#
VMR="VMRs/VMR_MSL39_v3.xlsx"
VMR_TSV="processed_accessions_b.tsv"
ACC_TSV="processed_accessions_b_acc_only.txt"
BATCH_N=1
OUT_DIR="refseq_${BATCH_N}"
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
## ENV
##
if [ -z "$(which conda 2>/dev/null)" ]; then
    echo "module load Anaconda3"
    module load Anaconda3
fi 

echo conda activate  conda/vmr_openpyxl3.working
conda activate  conda/vmr_openpyxl3.working

##
## PARSE MSL, extract genbank accession list
##
if [ $QUERY_ONLY -lt 1 ]; then
    
    if [[ $VMR -nt $VMR_TSV ]]; then 
	echo RUN ./VMR_to_fasta.py -verbose -VMR_file_name $VMR -mode VMR -ea B
	./VMR_to_fasta.py -verbose -VMR_file_name $VMR -mode VMR -ea B
    else
	echo "SKIP: $VMR_TSV up-to-date"
    fi
    wc -l $VMR_TSV

    # just accessions - remove ()'s
    cut -f 5   $VMR_TSV \
	| perl -pe 's/\([0-9.]+\)//g;' \
	> $ACC_TSV
    
    wc -l $ACC_TSV

    # split main file into small piecees
    mkdir -p $OUT_DIR
    split -l ${BATCH_N}  $ACC_TSV $PREFIX

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
