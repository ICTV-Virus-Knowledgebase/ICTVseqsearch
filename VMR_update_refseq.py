#!/usr/bin/env python3
#
# VMR_update_refseq.py
#
# Extract accessions and lineages from VMR.xls, query NCBI, build VMR.blast_db, query VMR.blast_db
#
# INPUT: 
#     -VMR_file_name VMRs/VMR_MSL39_v3.xlsx
#     -map           genbank_refseq_map_wide.txt
# ARGS: 
#     -      
# ITERMEDIATE FILES:
print("imports....")
import pandas as pd
import subprocess
import time
from urllib import error
import argparse
import numpy as np
import re
import sys
import os
import pathlib # for stem=basename(.txt)
# my utilities
from accession_utils import parse_seg_accession_list, merge_acc_dicts

parser = argparse.ArgumentParser(description="")

#setting arguments.
parser.add_argument('-verbose',help="printout details during run",action=argparse.BooleanOptionalAction)
parser.add_argument('-tmi',help="printout Too Much Information during run",action=argparse.BooleanOptionalAction)
parser.add_argument("-VMR_file_name",help="name of the VMR file to load.",default="VMRs/VMR_MSL39_v4/VMR_MSL39_v4c.xlsx")
parser.add_argument("-map_file_name",help="name TSV mapping genbank to refseq.",default="genbank_refseq_map_wide_fixed2.txt")
parser.add_argument("-ea",help="Fetch E or A records, or Both (Exemplars or AdditionalIsolates)", default="B")
args = parser.parse_args()

if args.ea.lower() != "b":
    print("Valid -EA mode not selected. Supported Options: B",file=sys.stderr)
    exit(1)
    
VMR_file_name_tsv = './vmr.tsv'
processed_accession_file_name ="./processed_accessions_"+args.ea.lower()+".tsv"


###############################################################################################################
# Loads excel from https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/ and puts it into a DataFrame
# NOTE: URL is incorrect. 
############################################################################################################### 
# DataFrame['column name'] = provides entire column
# DataFrame['column name'][0,1,2,3,4,5 etc] provides row for that column
# 
#
def load_VMR_data():
    if args.verbose: print("load_VMR_data()")
    if args.verbose: print("  opening", args.VMR_file_name)

    # Importing excel sheet as a DataFrame. Requires xlrd and openpyxl package
    try:
        # open excel file
        vmr_excel = pd.ExcelFile(args.VMR_file_name,engine='openpyxl')
        if args.verbose: print("\tOpened VMR Excel file: with {0} sheets: {1}".format(len(vmr_excel.sheet_names),args.VMR_file_name))

        # find first sheet matching "^VMR MSL"
        sheet_name = next((sheet for sheet in vmr_excel.sheet_names if re.match(r"^VMR MSL", sheet)), None)
        if args.verbose: print("\tFound sheet '{0}'.".format(sheet_name))

        if sheet_name is None:
            raise ValueError("No worksheet name matching the pattern '^VMR MSL' found.")
        else:
            raw_vmr_data = pd.read_excel(args.VMR_file_name,sheet_name=sheet_name,engine='openpyxl')
            if args.verbose: print("VMR data loaded: {0} rows, {1} columns.".format(*raw_vmr_data.shape))
            if args.verbose: print("\tcolumns: ",raw_vmr_data.columns)

            # list of the columns to extract from raw_vmr_data
            vmr_cols_needed = ['Isolate ID','Virus GENBANK accession','Virus REFSEQ accession']
            
            for col_name in list(raw_vmr_data.columns):
                if col_name in vmr_cols_needed:
                    print("    "+col_name+" [NEEDED]")
                else:
                    print("    "+col_name)
                    
            for col_name in vmr_cols_needed:
                if not col_name in list(raw_vmr_data.columns):
                    print("    "+col_name+" [!MISSING!]")

    except(FileNotFoundError):
        print("The VMR file specified does not exist! Make sure the path set by '-VMR_file_name' is correct.",file=sys.stderr)
    

    # save As TSV for diff'ing
    if os.path.exists(VMR_file_name_tsv) and os.path.getmtime(VMR_file_name_tsv) > os.path.getmtime(args.VMR_file_name):
        if args.verbose: print("  SKIP writing", VMR_file_name_tsv)
    else:
        if args.verbose: print("  writing", VMR_file_name_tsv)
        raw_vmr_data.to_csv(VMR_file_name_tsv,sep='\t', index=False)

    return raw_vmr_data

def parse_VMR_accessions(vmr_df):

    # Initialize an empty DataFrame to collect all the merged accessions
    vmr_accessions = pd.DataFrame()

    # Iterate over each row in the DataFrame
    for index, row in vmr_df.iterrows():
        isolate_id       = row['Isolate ID']
        genbank_acc_list = row['Virus GENBANK accession']
        refseq_acc_list  = row['Virus REFSEQ accession']

        # Parse the refseq and genbank accession lists
        refseq_accessions = parse_seg_accession_list(isolate_id, refseq_acc_list)
        genbank_accessions = parse_seg_accession_list(isolate_id, genbank_acc_list)

        # Merge the parsed dictionaries and get a DataFrame
        merged_accessions = merge_acc_dicts(isolate_id, genbank_accessions, refseq_accessions)

        # Append the resulting DataFrame to vmr_accessions
        vmr_accessions = pd.concat([vmr_accessions,merged_accessions], ignore_index=True)

    # write parsed accessions to TSV
    pd.DataFrame.to_csv(vmr_accessions, processed_accession_file_name, sep='\t')
    print("Wrote to ", processed_accession_file_name)
    if args.verbose: print("   Wrote {0} rows, {1} columns.".format(*vmr_accessions.shape))

    return vmr_accessions
    
    
###############################################################################################################
# Loads TSV RefSeq -> Genbank accession map
# Columsn:
############################################################################################################### 
# DataFrame['column name'] = provides entire column
# DataFrame['column name'][0,1,2,3,4,5 etc] provides row for that column
# 
#
def load_accession_map():
    if args.verbose: print("load_MAP_data()")
    if args.verbose: print("  opening", args.map_file_name)

    # Importing TSV as a DataFrame.
    try:
        # openfile
        raw_map_data = pd.read_csv(args.map_file_name,sep="\t")
        if args.verbose: print("MAP data loaded: {0} rows, {1} columns.".format(*raw_map_data.shape))
        if args.verbose: print("\tcolumns: ",raw_map_data.columns)

        # list of the columns to extract from raw_map_data
        cols_needed = ['#refseq','genbank','ncbi_taxid']

        for col_name in list(raw_map_data.columns):
            if col_name in cols_needed:
                print("    "+col_name+" [NEEDED]")
            else:
                print("    "+col_name)

        for col_name in cols_needed:
            if not col_name in list(raw_map_data.columns):
                print("    "+col_name+" [!MISSING!]")

    except(FileNotFoundError):
        print("The MAP ile specified does not exist! Make sure the path set by '-map_file_name' is correct.",file=sys.stderr)
    

    return raw_map_data


def main():
    if args.verbose: print("main()")

    print("# load VMR: ")
    vmr_df = load_VMR_data()

    print("# parse VMR accessions: ")
    accessions_df = parse_VMR_accessions(vmr_df)
    
    print("# load map: ")
    map_data = load_accession_map()

main()

if args.verbose: print("Done.")


