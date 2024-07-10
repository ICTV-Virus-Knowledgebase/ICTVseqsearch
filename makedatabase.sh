#!/usr/bin/env bash
#
# 20240117 runtime ~7h, RAM=402M
#
#SBATCH --job-name=ICTV_VMR_makeblastdb_e
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
#SBATCH --mem-per-cpu=30000
#
ACCESSION_TSV=processed_accessions_e.tsv
ALL_FASTA=./fasta_new_vmr/vmr_e.fa
SRC_DIR=$(dirname $ALL_FASTA)

echo "# concatenate individual fasta's into all.fa"
# if any fastas are updated, rebuild master fasta
if [[ ! -e "$ALL_FASTA" || "$(find $SRC_DIR -newer $ALL_FASTA|wc -l)" -gt 0 ]]; then
    echo "REBUILD $ALL_FASTA"
    rm $ALL_FASTA
    for FA in $(awk 'BEGIN{FS="\t";GENUS=5;ACC=3}(NR>1){print $GENUS"/"$ACC".fa"}' $ACCESSION_TSV); do
	echo "cat $SRC_DIR/$FA >> $ALL_FASTA"
	cat $SRC_DIR/$FA >> $ALL_FASTA
	wc -l $ALL_FASTA
    done
else
    echo "SKIP: $ALL_FASTA is up-to-date."
fi


echo "# Make the BLAST database"
if [ "$(which makeblastdb 2>/dev/null)" == "" ]; then 
    echo "module load BLAST"
    module load BLAST
fi

echo 'makeblastdb -in $ALL_FASTA -input_type "fasta" -title "ICTV VMR refseqs" -out "./blast/ICTV_VMR_e" -dbtype "nucl"'
makeblastdb -in $ALL_FASTA -input_type "fasta" -title "ICTV VMR refseqs" -out "./blast/ICTV_VMR_e" -dbtype "nucl"

echo "# Example usage:"
echo "# blastn -db ./blast/ICTV_VMR_e -query ./fasta_new_vmr/Eponavirus/MG711462.fa -out ./results/e/Eponavirus/MG711462.csv  -outfmt '7 delim=,'"
