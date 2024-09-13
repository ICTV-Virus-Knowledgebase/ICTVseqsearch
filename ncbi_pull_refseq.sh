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
#OUT_DIR="test_${BATCH_N}"
PREFIX="b_genbank.txt."
PREFIX_OUT="b_refseq.txt."
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

    # just accessions column
    # AND remove ()'s
    # AND remove version numebrs
    ./cutcols $VMR_TSV Accession \
	| perl -pe 's/\([0-9.]+\)//g;' \
	| perl -pe 's/\.[0-9]+$//g;' \
	> $ACC_TSV
    
    wc -l $ACC_TSV

    # split main file small files named by access number
    mkdir -p $OUT_DIR
    tail -n +2 $ACC_TSV | \
	awk -v OUT_DIR="$OUT_DIR" -v PREFIX="$PREFIX" '{subd = OUT_DIR"/"substr($1,length($1)-1,2); cmd="mkdir -p " subd; cmd  | getline; close(cmd); outf=subd"/"PREFIX""$1; print $1 > outf ; close(outf)}'

    # short test files
    head -2  $ACC_TSV  | grep -v Accession_IDs > $TEST1
    head -11  $ACC_TSV | grep -v Accession_IDs > $TEST10

fi

##
## QUERY NCBI
##

# ----------------------------------------------------------------------
# NCBI query for each file (TEST)
# --------------------------------------------------------------------
mkdir -p test
for file in $TEST0 $TEST1 $TEST10 ; do
    out=test/${file}.tsv
    err=test/${file}.err
    xml=test/${file}.xml
    echo "############### $file -> $out ##############"
    if [[ ! -e ${out} || ${file} -nt ${out} || $0 -nt ${out} ]]; then 
	echo "# QUERY NCBI: "
	#
	# check if we've fetched the record
	#
	if [[ ! -e ${xml} || ${file} -nt ${xml} ]]; then
	    echo "#      FETCH ${file}"
	    ( \
	      elink -input $file -db nucleotide -target nucleotide -name nuccore_nuccore_gbrs 2> $err \
		  || echo "ERROR $? : elink -input $file -db nucleotide -target nucleotide -name nuccore_nuccore_gbrs > $t1" > /dev/stderr \
	      ) \
		| \
		( \
		  efetch -format docsum > ${xml} \
		      || echo "ERROR $? : effect -format docsum > $xml"  > /dev/stderr \
		) 
	    # error check
	    grep -q "QUERY FAILURE" $err 2>/dev/null
	    if [ $? -eq 0 ]; then
		echo "ERROR: QUERY FAILURE: $file" 
		cat $err > /dev/stderr
	    fi
	else
	    echo "# SKIP FETCH ${file}"
	fi
	#
	# parse fetched data
	#
	if [[ ! -e ${out}|| ${xml} -nt ${out} || $0 -nt ${out} ]]; then
	    echo "#      PARSE ${file}"
	    ELEMENT_LIST="Caption AssemblyAcc Organism TaxId Slen MolType Topology Completeness Strand"
	    cat $xml | xtract -pattern DocumentSummary -sep "\t" -element $ELEMENT_LIST > $out \
		|| echo "ERROR $? : xtract -pattern DocumentSummary -sep '\t' -element $ELEMENT_LIST < $xml > $out">/dev/stderr
	    cat $out
	else
	    echo "# SKIP PARSE ${file}"
	fi
    else
	echo "# SKIP NCBI FETCH ${file}"
    fi
done

# ----------------------------------------------------------------------
# NCBI query for each file (REAL)
# ----------------------------------------------------------------------
for file in $OUT_DIR/*/${PREFIX}* ; do
    out=$(echo $file | perl -pe "s|$PREFIX|$PREFIX_OUT|;").tsv
    err=$(echo $file | perl -pe "s|$PREFIX|$PREFIX_OUT|;").err
    xml=$(echo $file | perl -pe "s|$PREFIX|$PREFIX_OUT|;").xml
    if [ $VERBOSE ]; then echo "## [$file => $out] ##  "; fi
    if [[ ! -e ${out} || ${file} -nt ${out} || $0 -nt ${out} ]]; then
	echo "# QUERY NCBI: "
	#
	# check if we've fetched the record
	#
	if [[ ! -e ${xml} || ${file} -nt ${xml} ]]; then
	    echo "#      FETCH ${file}"
	    ( \
	      elink -input $file -db nucleotide -target nucleotide -name nuccore_nuccore_gbrs 2> $err \
		  || echo "ERROR $? : elink -input $file -db nucleotide -target nucleotide -name nuccore_nuccore_gbrs > $t1" > /dev/stderr \
	      ) \
		| \
		( \
		  efetch -format docsum > ${xml} \
		      || echo "ERROR $? : effect -format docsum > $xml"  > /dev/stderr \
		) 
	    # error check
	    grep -q "QUERY FAILURE" $err 2>/dev/null
	    if [ $? -eq 0 ]; then
		echo "ERROR: QUERY FAILURE: $file" 
		cat $err > /dev/stderr
	    fi
	else
	    echo "# SKIP FETCH ${file}"
	fi
	#
	# parse fetched data
	#
	if [[ ! -e ${out}|| ${xml} -nt ${out} || $0 -nt ${out} ]]; then
	    echo "#      PARSE ${file}"
	    ELEMENT_LIST="Caption AssemblyAcc Organism TaxId Slen MolType Topology Completeness Strand"
	    cat $xml | xtract -pattern DocumentSummary -sep "\t" -element $ELEMENT_LIST > $out \
		|| echo "ERROR $? : xtract -pattern DocumentSummary -sep '\t' -element $ELEMENT_LIST < $xml > $out">/dev/stderr
	    cat $out
	else
	    echo "# SKIP PARSE ${file}"
	fi
    else
	echo "# SKIP NCBI FETCH ${file}"
    fi
done

echo "#"
echo "# merge output files"
echo "#"
#sort -k2,2 $OUT_DIR/*/${PREFIX_OUT}* > $OUT_MAP # too many files
echo "#" find $OUT_DIR -name "${PREFIX_OUT}*.tsv" -exec cat {} + | sort -k2,2 > $OUT_MAP
find $OUT_DIR -name "${PREFIX_OUT}*.tsv" -exec cat {} + | sort -k2,2 > $OUT_MAP

echo "#"
echo "# scan for errors"
echo "#"
echo ./list_errors.sh $OUT_DIR \"\|\"
./list_errors.sh $OUT_DIR \|

echo "QC $VMR"
echo "   GENBANK: $(wc -l $VMR_TSV)"
echo "   REFSEQ:  $(wc -l $OUT_MAP)"
echo "   ERRORS:  $(wc -l $OUT_DIR/error_list.txt)"


echo "#"
echo "# merge RefSeq into VMR export"
echo "#"
echo ./update_vrm_refseq.py                                                   
./update_vrm_refseq.py                                                   

