#!/usr/bin/env bash
#
# test_e_accessions.sh [-g genus] [-m mode]
#
# mode=default,blastn,blastn16,blastn20,blastn24,blastn-short
#
# query each "A" (additional) sequence against the blastdb of "E" sequences. 
# tabulate matches and mis-matches
#
# (This would be much faster if we merged all the query fastas)
#
#
# 20240122 runtime 
#
#SBATCH --job-name=ICTV_BLAST_A_fastas_vs_E_db
#SBATCH --output=logs/log.%J.%x.out
#SBATCH --error=logs/log.%J.%x.err
#
# Number of tasks needed for this job. Generally, used with MPI jobs
# Time format = HH:MM:SS, DD-HH:MM:SS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000  # in Megabytes
# 
# partition and time
#SBATCH --partition=amd-hdr100 --time=00-12:00:00
##SBATCH --partition=amd-hdr100 --time=06-06:00:00
##SBATCH --partition=medium --time=40:00:00
#

QUERY_DIR=./fasta_new_vmr_a
GENUS="*"
if [[ "$1" == "-g" && ! -z "$2" ]]; then GENUS=$2; shift 2; fi
echo "GENUS:       $GENUS"
BLAST_DB=./blast/ICTV_VMR_e
echo "BLAST_DB:    $BLAST_DB"
BLAST_OUT_FMT="-outfmt '7 delim=,'"; BLAST_OUT_SUFFIX="raw.txt"  # scot original
BLAST_OUT_FMT="-outfmt 5 -max_target_seqs 10 -max_hsps 10"; BLAST_OUT_SUFFIX="hit10hsp10.xml"
BLAST_OUT_FMT="-outfmt 5"; BLAST_OUT_SUFFIX="full.xml"
BLAST_OUT_FMT="-outfmt 5 -max_target_seqs 10"; BLAST_OUT_SUFFIX="hit10.xml"
BLAST_OUT_FMT="-outfmt 5 -max_target_seqs 20"; BLAST_OUT_SUFFIX="hit20.xml"
BLAST_OUT_FMT="-outfmt 5 -max_target_seqs 60"; BLAST_OUT_SUFFIX="hit60.xml"
# 
# BLAST_ARGS
# 
# The default blastn command runs megablastn with a default word size of 28
# '-task blastn'       uses a default word size of 11.
# '-task blastn-short' uses a default word size of 7
BLAST_ABBREV="default"
#BLAST_ABBREV="blastn"
#BLAST_ABBREV="short"
if [[ "$1" == "-m" && ! -z "$2"   ]]; then 
    BLAST_ABBREV=$2
    echo BLAST_ABBREV=$BLAST_ABBREV
    shift 2
fi
if [[ $BLAST_ABBREV == "default"  ]]; then BLAST_ARGS=""; fi
if [[ $BLAST_ABBREV == "blastn"   ]]; then BLAST_ARGS="-task blastn"; 
else 
    if [[ "$BLAST_ABBREV" == blastn* ]]; then 
	BLAST_ARGS="-task blastn -word_size $(echo $BLAST_ABBREV|sed 's/blastn//')"
	echo "BLAST_ARGS=\"$BLAST_ARGS\""
    fi
fi
if [[ $BLAST_ABBREV == "short"    ]]; then BLAST_ARGS="-task blastn-short"; fi

RESULTS_DIR=./results/$BLAST_ABBREV/a
echo "RESULTS_DIR: $RESULTS_DIR"
mkdir -p $RESULTS_DIR

if [ -z "$(which blastn 2>/dev/null)" ]; then
    echo "#module load BLAST"
    echo module load BLAST
    module load BLAST
    echo "which blastn: " $(which blastn)
fi


ACCESSION_TSV=./processed_accessions_a.tsv
echo "#"
echo "# enumerate A seqs"
echo "# "
echo "# scan ACCESSION_TSV=$ACCESSION_TSV"
echo "# "
A_FASTAS=$(awk -v DIR=$QUERY_DIR -v TARGET_GENUS="$GENUS" 'BEGIN{FS="\t";GENUS=5;ACC=3}(NR>1&&(TARGET_GENUS=="*"||TARGET_GENUS==$GENUS)){print DIR "/" $GENUS "/" $ACC ".fa"}' $ACCESSION_TSV)
#A_FASTAS=$(find $QUERY_DIR -name "*.fa")
echo '# found ' $(wc -l $ACCESSION_TSV) " A accessions"
if [ "$GENUS" != "*" ]; then 
    echo '# FILTERD TO ' $(echo $A_FASTAS | sed -e 's/ /\n/g'| wc -l ) " GENUS=$GENUS A accessions"
fi

echo "#"
echo "# run blast"
echo "#" 
for QUERY in $A_FASTAS; do
    
    GENUS=$(basename $(dirname $QUERY))
    ACCESS=$(basename $QUERY .fa)
    RESULT_DIR=$RESULTS_DIR/$GENUS
    RESULTS_FILE_RAW="$RESULT_DIR/$ACCESS.$BLAST_OUT_SUFFIX"
    RESULTS_FILE_SORT_BIT="$RESULT_DIR/$ACCESS.bitscore.csv"
    RESULTS_FILE_SORT_EV="$RESULT_DIR/$ACCESS.evalue.csv"

    # 
    # run BLAST
    #
    if [[ -e $RESULTS_FILE_RAW && $RESULTS_FILE_RAW -nt $QUERY ]]; then
	echo "${GENUS}/${ACCESS}: SKIP BLAST"
    else
	echo -n "${GENUS}/${ACCESS}: blast " 
	mkdir -p $RESULT_DIR
	blastn $BLAST_ARGS -db $BLAST_DB -query $QUERY -out ${RESULTS_FILE_RAW} ${BLAST_OUT_FMT}   
	echo -n $(grep "hits found" $RESULTS_FILE_RAW)
	echo " # blastn $BLAST_ARGS -db $BLAST_DB -query $QUERY -out ${RESULTS_FILE_RAW} ${BLAST_OUT_FMT}"
    fi

    # #
    # # SORT blast output: BITSCORE
    # # 
    # if [[ -e $RESULTS_FILE_SORT_BIT && $RESULTS_FILE_SORT_BIT -nt $RESULTS_FILE_RAW ]]; then
    # 	echo "${GENUS}/${ACCESS}: SKIP SORT bitscore"
    # else
    # 	echo "${GENUS}/${ACCESS}: sort by bitscore " 
    # 	(head -n 5 $RESULTS_FILE_RAW; \
    #      grep -v "^#" $RESULTS_FILE_RAW | sort -t , -k12,12nr; \
    # 	 tail -n 1 $RESULTS_FILE_RAW ) \
    #      > $RESULTS_FILE_SORT_BIT
    # fi

    # #
    # # SORT blast output: EVALUE
    # # 
    # if [[ -e $RESULTS_FILE_SORT_EV && $RESULTS_FILE_SORT_EV -nt $RESULTS_FILE_RAW ]]; then
    # 	echo "${GENUS}/${ACCESS}: SKIP SORT evalue"
    # else
    # 	echo "${GENUS}/${ACCESS}: sort by evalue " 
    # 	(head -n 5 $RESULTS_FILE_RAW; \
    #      grep -v "^#" $RESULTS_FILE_RAW | sort -t , -k11,11n -k12,12rn; \
    # 	 tail -n 1 $RESULTS_FILE_RAW ) \
    #      > $RESULTS_FILE_SORT_EV
    # fi

done
