#!/usr/bin/env bash
#
#  download latest VMR xlsx's from box
#  parse latest VMR (or $1)

# ---------------------------------------------------------------------- 
# environment
# ---------------------------------------------------------------------- 
#Set your environment here
if [ -z "$(which conda 2>/dev/null)" ]; then
    echo module load Anaconda3
    module load Anaconda3
fi
if [[ "$(which python 2>/dev/null)" != *$PWD/conda* ]]; then
    echo conda activate conda/vmr_openpyxl3
    conda activate conda/vmr_openpyxl3
fi

echo "#"
echo "# download VMRs from box into ./VMRs/"
echo "#"
echo  "#./pull_VMRs_from_box.sh"
$SRUN ./pull_VMRs_from_box.sh

echo "#"
echo "# get most recent MSL.xlsx file"
echo "#"

VMR_XLSX=$(find VMRs -name "VMR_MSL*.xlsx" -exec ls -rt {} +| tail -1)
if [ ! -z "$1" ]; then 
    VMR_XLSX=$1
fi
echo VMR_XLSX=$VMR_XLSX

echo "#"
echo "# parse VMR "
echo "#"
echo  ./VMR_to_fasta.py -verbose -mode VMR -VMR_file_name $VMR_XLSX
$SRUN ./VMR_to_fasta.py -verbose -mode VMR -VMR_file_name $VMR_XLSX

