#!/usr/bin/env bash
#
# Summarize status of an NCBI/RefSeq query/cache dir
#
if [[ "$1" == -h* ]]; then 
    echo "# Syntax: ./status.sh [dir]"
    exit 0
fi
#
# Parse args
#
DIR="." ; if [ ! -z "$1" ]; then DIR="$1"; shift;fi

#
# tabulate query/parse status
#
(
    echo "dir : gb xml tsv err"
    for d in $(find $DIR -name "??" -type d -maxdepth 1|sort); do
	echo $(basename $d) : \
	     $(ls $d/* 2>/dev/null | egrep -v "(xml|tsv)$" | wc -l) \
	     $(ls $d/*.xml 2>/dev/null | wc -l) \
	     $(ls $d/*.tsv 2> /dev/null| wc -l ) \
	     $(find $d -name "b_refseq*.err" -size +0 | wc -l )
    done
) 2> /dev/null
