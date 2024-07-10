#!/usr/bin/env bash
#
# query each "A" (additional) sequence against the blastdb of "E" sequences. 
#
# split by genus
#
#
MODES="blastn10 blastn11 blastn13 blastn16 blastn20 blastn24 blastn13 blastn15   blastn16 blastn17  blastn18 blastn22 blastn28 blastn9 "
MODES="blastn10 blastn11 blastn13"
MODES="blastn7"
MODES="blastn10 blastn7"
MODES="blastn10"
if [ ! -z "$*" ]; then
    MODES="$*"
fi
echo MODES=$MODES

GENERA=$(cd fasta_new_vmr_a/;ls -ls|grep drwx|sed 's/.* //g')
echo GENERA N=$(echo $GENERA|wc -w)

for MODE in $MODES; do
    echo "MODE=$MODE"
    for GENUS in $GENERA; do
	echo -n "   [$MODE] [$GENUS] "
	sbatch --job-name ICTV_BLAST_A_fastas_vs_E_db_${MODE}_$GENUS ./test_e_accessions.sh -g $GENUS -m $MODE
	# prevent SLURM overload, which looks like
	# sbatch: error: Batch job submission failed: Socket timed out on send/recv operation
	sleep 2
    done
done

