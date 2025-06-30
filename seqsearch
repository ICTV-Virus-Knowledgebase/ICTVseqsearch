#!/usr/bin/env python3
#
# taxonomically classify a novel nucleotide sequence
#

import argparse
import csv
import re
import pandas as pd
import json
import os
import sys
import subprocess
import glob
import numpy as np
#
# parse ARGS
#
parser = argparse.ArgumentParser(description="")
parser.add_argument('-verbose',help="print details during run",action=argparse.BooleanOptionalAction)
parser.add_argument('-indir',help="directory for fasta files with NUCLEOTIDE seqeunces to classify",default="seq_in")
parser.add_argument('-outdir',help="output directory for tax_results.json",default="tax_out")
parser.add_argument('-json',help="output json filename",default="tax_results.json")
parser.add_argument('-html',help="output  html base filename",default="tax_results")
parser.add_argument('-blastdb',help="input blast base filename",default="./blast/ICTV_VMR_b")
parser.add_argument('-db', help="BLAST database file", default="blast/ICTV_VMR_b")
parser.add_argument('-info', action='store_true', help="Display BLAST database info")

args = parser.parse_args()

print("Hello World")
print("IN_DIR: ", args.indir)
print("OUT_DIR:", args.outdir)
print("OUT_JSON:", args.json)
print("OUT_HTML:", args.html)

#
# read version number
#
with open("version_git.txt") as verfile:
    git_version = verfile.readline().strip()
print("VERSION:", git_version )

#System command to get blast db info "blastdbcmd -db blast/ICTV_VMR_b -info"
db_title_command = ["blastdbcmd", "-db", args.db, "-info"]
db_title_output = subprocess.run(db_title_command, capture_output=True, text=True, check=True)

#Split the output into lines and take the first 4 lines with splicing
db_title_len= db_title_output.stdout.splitlines()
first_4_lines = db_title_len[:4]
# Joins the first 4 lines by the tab characters
if first_4_lines:
    first_4_lines= "\t".join(first_4_lines)
    #Removes any tabs in the first 4 lines and replaces them with a single space
    remove_tabs = re.sub(r'\t+', ' ', first_4_lines)
    first_4_lines = remove_tabs.strip()
    remove_words= first_4_lines[9:]  # remove the first 9 characters "Database: "
    first_4_lines = remove_words.strip()





json_hits = {
    "program_name": os.path.basename(sys.argv[0]),
    "version"  : git_version,
    "database_name": args.blastdb,
    "database_title": first_4_lines,
    "input_dir": args.indir,
    "errors": "",
    "results":[]          
}

tsv_hits = None

rankmap = {
    "Realm": "realm",
    "Family": "family",
    "Subfamily": "subfamily",
    "Phylum": "phylum",
    "Class": "class",
    "Order": "order",
    "Genus": "genus",
    "Species": "species",
    "Kingdom": "kingdom",
    "Subkingdom": "subkingdom",
    "Subphylum": "subphylum",
    "Subrealm": "subrealm",
    "Subclass": "subclass",
    "Suborder": "suborder",
    "Subgenus": "subgenus",
}


GroupedHits = {}






# for eden: get list of fasta files in input directory (.fa/.fna.fasta) 
if args.verbose:
    print("Files in input directory: ", os.listdir(args.indir))
# get list of fasta files
fasta_files = []
for  infilename in os.listdir(args.indir):
    if infilename.endswith(".fa") or infilename.endswith(".fna") or infilename.endswith(".fasta"):
            fasta_files.append (infilename)
print("FASTA files found: ", fasta_files)

 #parsed_data_df is df for the processed_accessions_b.tsv file
parsed_data_df= pd.read_csv("processed_accessions_b.fa_names.tsv", sep="\t", header=0)


# iterate over fasta files
for infilename in sorted(fasta_files):
    # read in fasta file
    print("\n\n----------------------------------------------------------------------")
    print("Processing fasta file: ", infilename)
    fasta_file_path = os.path.join(args.indir, infilename)
    #create output file paths
    blastasn_output_filepath= os.path.join(args.outdir, infilename + ".asn")
    blastcsv_output_filepath= os.path.join(args.outdir, infilename + ".csv")
    blasthtml_output_filepath= os.path.join(args.outdir, infilename + ".html")

    #blast output format commands.

    blastasn_cmd = ["blastn", "-db", args.blastdb , "-query", fasta_file_path , "-out", blastasn_output_filepath, "-outfmt", "11"]
    blastcsv_cmd = ["blast_formatter", "-archive", blastasn_output_filepath , "-out", blastcsv_output_filepath, "-outfmt", "10"]
    blasthtml_cmd = ["blast_formatter", "-archive", blastasn_output_filepath,  "-out", blasthtml_output_filepath, "-html" ]
    print("ASN Format Command:", " ".join(blastasn_cmd), "\n", "CSV Format Command:", " ".join(blastcsv_cmd), "\n", "HTML Format Command:", " ".join(blasthtml_cmd))

    db_output_asn = subprocess.run(blastasn_cmd, capture_output=True, text=True)
    print("Blast .asn output files:", blastasn_output_filepath)
    db_output_csv = subprocess.run(blastcsv_cmd, capture_output=True, text=True)
    print("Blast .csv output files:", blastcsv_output_filepath)
    db_output_html = subprocess.run(blasthtml_cmd, capture_output=True, text=True)
    print("Blast .html output files:", blasthtml_output_filepath)




    
    # read in csv file and set headers to extract data
    
    try:
        raw_data = pd.read_csv(blastcsv_output_filepath)
        # Process the file only if it's successfully read
        #df is df for the blast csv file
        

        df = pd.read_csv(blastcsv_output_filepath, header=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
        print( "Loading output file", blastcsv_output_filepath)
        print( "Number of rows in blast csv file: ", len(df))
        if args.verbose: print(df)
        
        print("Number of rows in processed_accessions_b.tsv file: ", len(parsed_data_df))
        #initializes the sBaseAccession column in df that is used for joining
        df["sBaseAccession"]= df["sseqid"].astype(str).str.replace(r'^.*-(\w[^.]*)\..*$', r'\1',regex=True)
        print(df["sBaseAccession"])
        #merging the two dataframes by the sBaseAccession column and the Accession column
        merge_dfs= pd.merge(df,parsed_data_df, left_on="sBaseAccession", right_on= "Accession", how="left")
        if args.verbose: print(merge_dfs)

        # check if some blast results didn't match to known accessions
        merge_missing = merge_dfs["Species"].isna().sum()
        if merge_missing > 0: 
            print(f"INTERNAL ERROR: {merge_missing} out of {len(df)} sseqid.base_accession's didn't match to processed_accessions table.")
            print("INTERNAL ERROR: check if your blastdb and processed_accesions table are out of sync.")

        if len(df) < len(merge_dfs):
            print("INTERNAL ERROR: merging blast results to processed_accessions added extra rows.")
            print("INTERNAL ERROR: check if your blastdb and processed_accesions table are out of sync.")
            sys.exit(1)

        json_hits["results"]= list(GroupedHits.values())
    
                    
        #iterate over the merged dataframe to extract the data to hit dict
        for index, row in merge_dfs.iterrows():
            
        
            hit= {}           # append to hit
            hit["input_seq"] = infilename.split(".")[0]
            hit["evalue"]= row["evalue"]
            hit["bitscore"]= row["bitscore"]
            hit["qseqid"]= row["qseqid"]
            hit["sseqid"]= row["sseqid"]
            hit["ICTV_ID"]= row["ICTV_ID"]
            hit["isolate_ID"]= row["Isolate_ID"]
            hit["exemplar_additional"]= row["Exemplar_Additional"]
            hit["virus_names"]= row["Virus_Names"]



            #Turns empty NaN values to "Null"
            if isinstance(row["Start_Loc"], float) and np.isnan(row["Start_Loc"]):
                hit["start_loc"] = None
            else:
                hit["start_loc"] = int(row["Start_Loc"])
            if isinstance(row["End_Loc"], float) and np.isnan(row["End_Loc"]):
                hit["end_loc"] = None
            else:
                hit["end_loc"] = int(row["End_Loc"])
    
            #establishing dict in hit dict
            hit["sseqid_lineage"]= {}
            for src_key, target_key in rankmap.items():
                hit["sseqid_lineage"][target_key] = row[src_key]

                if isinstance(hit["sseqid_lineage"][target_key], float) and np.isnan(hit["sseqid_lineage"].get(target_key)):
                    hit["sseqid_lineage"][target_key] = None

            

            #Turns empty NaN values in the Segment_Name column into "Null"
            if isinstance(row["Segment_Name"], float) and np.isnan(row["Segment_Name"]):
                hit["segmentname"] = None
            else:
                hit["segmentname"] = row["Segment_Name"]
            #splitting the sseqid row into two parts
            if row["sseqid"]:
                hit["sseqid_accession"] = row["sseqid"].replace(row["sseqid"].rsplit("-", 1)[0] + "-", "")

           

            #splitting the qseqid row into two parts to just show the qseqid without everything after the first #
            delimiter= r"[#]"
            hit["qseqid"]= re.split(delimiter, hit["qseqid"])[0]

            if blastcsv_output_filepath not in GroupedHits:
                GroupedHits[blastcsv_output_filepath] = {
                    "input_file": infilename,
                    "blast_csv": re.sub(r".*/", "", blastcsv_output_filepath),  # just the filename
                    "blast_html": re.sub(r".*/", "", blasthtml_output_filepath),
                    "status": "HITS",
                    "errors": "",
                

                    "hits": []
                }
            GroupedHits[blastcsv_output_filepath]["hits"].append(hit)
            json_hits["results"]= list(GroupedHits.values())

            while json_hits["results"] is False:


            

            
                del hit["sseqid_lineage"]

                for src_key, target_key in rankmap.items():
                    hit[target_key] = row[src_key]
                    if isinstance(hit[target_key], float) and np.isnan(hit[target_key]):
                        hit[target_key] = None




# Write TSV file
                pd.DataFrame(GroupedHits[blastcsv_output_filepath]["hits"]).to_csv(
                    os.path.join(args.outdir, f"{os.path.basename(blastcsv_output_filepath)}_hits.tsv"))


                
                
                
                
                print("Wrote tsv file: ", os.path.join(args.outdir, blastcsv_output_filepath + "_hits.tsv"))


           
        
                    



            
            
            
    except pd.errors.EmptyDataError:
        print(blastcsv_output_filepath, "is empty and has been skipped.", "status: NO_HITS")







        


            
            
            
    
    
            





# Ensure the output directory exists
os.makedirs(args.outdir, exist_ok=True)

# json job output/summary
json_outpath=os.path.join(args.outdir, args.json)
with open(json_outpath, "w") as outfile:
    json.dump(obj=json_hits, fp=outfile, indent=4)
print("Wrote: ", json_outpath)
outfile.close()



exit(0)


