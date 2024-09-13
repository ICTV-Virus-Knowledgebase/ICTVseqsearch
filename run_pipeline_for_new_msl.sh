#!/usr/bin/env bash
#
# given a new MSL
#  download from box
#  parse
#  pull fasta from NCBI & reformat
#  build E record blastdb
#  blast A records against E-blastdb
#  tabulate results
#
#SBATCH --job-name=ICTV_VMR_pipeline
#SBATCH --output=logs/log.%J.%x.out
#SBATCH --error=logs/log.%J.%x.err
#
# Number of tasks needed for this job. Generally, used with MPI jobs
# Time format = HH:MM:SS, DD-HH:MM:SS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000  # in Megabytes
# 
# partition and time
#SBATCH --partition=amd-hdr100 --time=06-06:00:00
##SBATCH --partition=amd-hdr100 --time=00-12:00:00
##SBATCH --partition=medium --time=40:00:00

# ---------------------------------------------------------------------- 
# default settings / arg parsing
# ---------------------------------------------------------------------- 
BLAST_MODES="blastn10 blastn11 blastn13"
BLAST_MODES="blastn10"
echo BLAST_MODES=$BLAST_MODES

# ---------------------------------------------------------------------- 
# better step logging if sbatch'ed
# ---------------------------------------------------------------------- 
SRUN=""
if [ ! -z "$SLURM_JOB_ID" ]; then 
    SRUN=srun
    echo "SLURM_JOB_ID=$SLURM_JOB_ID"
    echo "seff $SLURM_JOB_ID"
    echo "sacctx -j $SLURM_JOB_ID"
fi
#SRUN="echo DEBUGGING: "
echo "SRUN='$SRUN'"

# ---------------------------------------------------------------------- 
# environment
# ---------------------------------------------------------------------- 
#Set your environment here
if [ -z "$(which conda 2>/dev/null)" ]; then
    echo module load Anaconda3
    module load Anaconda3
fi
if [[ "$(which python 2>/dev/null)" != *$PWD/conda* ]]; then
    if [[ -e "conda/vmr_openpyxl3.working" ]]; then
        echo conda activate conda/vmr_openpyxl3_working
        conda activate conda/vmr_openpyxl3_working
    else
        echo conda activate conda/vmr_openpyxl3
        conda activate conda/vmr_openpyxl3
    fi
fi

# ---------------------------------------------------------------------- 
# start work
# ---------------------------------------------------------------------- 
echo "#"
echo "# download VMRs from box into ./VMRs/"
echo "#"
echo  "#./pull_VMRs_from_box.sh"
$SRUN ./pull_VMRs_from_box.sh

echo "#"
echo "# get most recent MSL.xlsx file"
echo "#"

VMR_XLSX=$(ls -rt VMRs/VMR_MSL*.xlsx | tail -1)
if [ ! -z "$1" ]; then 
    VMR_XLSX=$1
fi
echo VMR_XLSX=$VMR_XLSX

echo "#"
echo "# parse VMR "
echo "#"
echo  ./VMR_to_fasta.py -verbose -mode VMR -ea b -VMR_file_name $VMR_XLSX
$SRUN ./VMR_to_fasta.py -verbose -mode VMR -ea b -VMR_file_name $VMR_XLSX

echo "#"
echo "# pull and format fastas"
echo "#"
echo "# remove formated fastas"
echo '#SKIP: find . -name "*.fa" -exec rm {} +'
#$SRUN find . -name "*.fa" -exec rm {} +

echo "#SKIP: downlaod and format E fastas - in serial"
#$SRUN ./download_fasta_files $VMR_XLSX

echo "#SKIP: downlaod and format A fastas - in serial"
#$SRUN ./download_fasta_files_a  $VMR_XLSX

echo "#"
echo "# build blast db: 3-4 minutes"
echo "#"
echo "makedatabase.sh"
$SRUN makedatabase.sh

echo "#"
echo "# test blast A records vs E_blastdb: 30-40 minutes, in parallel"
echo "#"
echo "sbatch_modes_by_genus_test_e.sh $BLAST_MODES"
$SRUN sbatch_modes_by_genus_test_e.sh $BLAST_MODES

echo "# wait for parallel blast jobs to complete"
echo "squeue_email ICTV_BLAST_A_fastas_vs_E_db_ done for $0"
$SRUN squeue_email ICTV_BLAST_A_fastas_vs_E_db_ done for $0

for BLAST_MODE in $BLAST_MODES; do

    echo "#"
    echo "# tabulate results"
    echo "#"
    echo  ./tabulate_test_results.py --blast $BLAST_MODE
    $SRUN ./tabulate_test_results.py --blast $BLAST_MODE

    echo "#"
    echo "# summarize results"
    echo "#"
    echo "cd results/"
    pushd results/
    echo "./extract_results.sh -m $BLAST_MODE -s bitscore"
    $SRUN ./extract_results.sh -m $BLAST_MODE -s bitscore
    popd
done

echo "#" 
echo "# announce done"
echo "#"
(
    echo BLAST_MODES=$BLAST_MODES
    echo "SLURM_JOB_ID=$SLURM_JOB_ID"
    echo "sacctx -j $SLURM_JOB_ID"
    echo "run in $PWD"
    echo " "
    echo "SUMMARY (table):"
    (for BLAST_MODE in $BLAST_MODES; do echo -e "$VMR_XLSX,$BLAST_MODE,$(cat ./$BLAST_MODE/a/summary.bitscore.mismatch_totals.tsv)"; done) | column -t -s , 
    echo " "
    echo "SUMMARY (CSV):"
    (for BLAST_MODE in $BLAST_MODES; do echo -e "$VMR_XLSX,$BLAST_MODE,$(cat ./$BLAST_MODE/a/summary.bitscore.mismatch_totals.tsv)"; done ) 
) \
| tee ../run_pipeline_for_new_msl.email.txt \
| $SRUN mail -s "$0 $* COMPLETED" $USER@uab.edu


