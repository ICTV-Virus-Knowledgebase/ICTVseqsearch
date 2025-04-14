#!/usr/bin/env bash
#
# 20240117 runtime ~7h, RAM=402M
#
#SBATCH --job-name=ICTV_VMR_pull_fastas_from_ncbi
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
#SBATCH --mem-per-cpu=3
#
EA=b
ACCESSION_TSV=processed_accessions_$EA.tsv
ALL_FASTA=./fasta_new_vmr_$EA/vmr_$EA.fa
SRC_DIR=$(dirname $ALL_FASTA)
BLASTDB=./blast/ICTV_VMR_$EA
FIRST_FASTA=$(awk 'BEGIN{FS="\t";GENUS=7;ACC=5}(NR>1){print $GENUS"/"$ACC".fa"}' $ACCESSION_TSV|head -1)
OUT_FILEPATH=$(awk 'BEGIN{FS="\t";GENUS=7;ACC=5}(NR>1){print $GENUS"/"$ACC}' $ACCESSION_TSV|head -1)

echo EA=$EA
echo ACCESSION_TSV=$ACCESSION_TSV
echo ALL_FASTA=$ALL_FASTA
echo SRC_DIR=$SRC_DIR
echo BLASTDB=$BLASTDB

if [[ -z "$(which conda 2>/dev/null)" ]]; then
    echo module load Anaconda3
    module load Anaconda3
fi

conda activate conda/vmr_openpyxl3

#
# download fastas
#
./VMR_to_fasta.py -mode fasta -ea b -email $USER@uab.edu -v

