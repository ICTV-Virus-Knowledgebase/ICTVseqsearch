#!/usr/bin/env python3
#
# Generate VMR RefSeq strings, and SQL updates to load them in the db
#
# INPUT: processed_accessions_e.tsv
#    Species Accession_IDs segment Sort Original_Accession_String
# INPUT: genbank_refseq_map.txt
#    REFSEQ GENBANK NCBI_ISOLATE_NAME    
#
# Update will use Original_Accession_String as the PK into VMR table.
# It will update both the Genbank and RefSeq accession strings
#
import argparse
import pandas
from pprint import pprint

parser = argparse.ArgumentParser(description="")

#setting arguments.
parser.add_argument('-verbose',help="printout details during run",action=argparse.BooleanOptionalAction, default=True)
parser.add_argument("-vmr_accessions",help="TSV: VMR accessions, one per line, with segment name", default="processed_accessions_b.tsv")
parser.add_argument("-refseq_map",help="TSV: maps Genbank accessios to matching RefSeq accessions", default="genbank_refseq_map.txt")
parser.add_argument("-sql",help="output filename for update SQL.",default="vmr_genbank_refseq_update.sql")
parser.add_argument("-tsv",help="output filename for updated refseq TSV", default="vmr_genbank_refseq_update.tsv")
args = parser.parse_args()

print("###############################################################################################################")
print("# Loads excel inputs")
print("###############################################################################################################")

if args.verbose: print("Loading: ", args.vmr_accessions)
vmr_proc_df = pandas.read_csv(args.vmr_accessions, sep="\t")
if args.verbose: print("   Read {0} rows, {1} columns.".format(*vmr_proc_df.shape))
if args.verbose: print("   columns:", list(vmr_proc_df.columns))

if args.verbose: print("Loading: ", args.refseq_map)
# name "refseq" Accession_IDs, so column names match for merge
#refseq_colnames=['refseq', 'Accession_IDs', 'ncbi_virus_name']
refseq_colnames=['refseq', 'genbank', 'ncbi_virus_name']
refseq_df = pandas.read_csv(args.refseq_map, sep="\t", header=None, names=refseq_colnames)
if args.verbose: print("   Read {0} rows, {1} columns.".format(*refseq_df.shape))
if args.verbose: print("   columns:", list(refseq_df.columns))

print("###############################################################################################################")
print("# clean data ")
print("###############################################################################################################")

#
# remove .# suffix from accession numbers
#
vmr_proc_df["genbank"] = vmr_proc_df.Accession_IDs.replace(to_replace='\\.[-0-9]+$', regex=True, value='')
refseq_df[  "genbank"] =   refseq_df.genbank.replace(      to_replace='\\.[-0-9]+$', regex=True, value='')

print("###############################################################################################################")
print("# merge data frames")
print("###############################################################################################################")

vmr_refseq_merge=pandas.merge(vmr_proc_df, refseq_df, on=['genbank'], how="outer", indicator=True)

vmr_refseq_df=vmr_refseq_merge.query('_merge!="right_only"')
bad_refseq_df=vmr_refseq_merge.query('_merge=="right_only"')

bad_fname = "vmr_genbank_refseq_unmatched.tsv"
if args.verbose: print("Writing: ", bad_fname)
pandas.DataFrame.to_csv(bad_refseq_df, bad_fname, sep="\t")
if args.verbose: print("   Read {0} rows, {1} columns.".format(*bad_refseq_df.shape))

if args.verbose: print("Writing: ", args.tsv)
pandas.DataFrame.to_csv(vmr_refseq_df, args.tsv, sep="\t")
if args.verbose: print("   Read {0} rows, {1} columns.".format(*vmr_refseq_df.shape))

print("###############################################################################################################")
print("# tabulate segments")
print("###############################################################################################################")


# Define a function to concatenate segment and refseq values
def format_seg_acc(seg,acc):
    if seg != seg: 
        return str(acc)
    return str(seg)+': '+str(acc)

# filter and format
x = vmr_refseq_df.dropna(subset=['refseq'])
x['seg_refseq']= x.apply(lambda r: format_seg_acc(r['segment'],r['refseq']), axis=1)
x['seg_genbank']=x.apply(lambda r: format_seg_acc(r['segment'],r['genbank']),axis=1)

# Group by Original_Accession_String and apply the concatenation function
groupby_columns = ['Species','Exemplar_Additional','Original_Accession_String']
x_new_refseq=    x[[*groupby_columns,'seg_refseq' ]].groupby(by=groupby_columns)['seg_refseq' ].apply(lambda v: '; '.join(v)).reset_index()
x_new_refseq_ct= x[[*groupby_columns,'seg_refseq' ]].groupby(by=groupby_columns)['seg_refseq' ].apply(lambda v: len(v))
x_new_genbank=   x[[*groupby_columns,'seg_genbank']].groupby(by=groupby_columns)['seg_genbank'].apply(lambda v: '; '.join(v)).reset_index()
x_new_genbank_ct=x[[*groupby_columns,'seg_genbank']].groupby(by=groupby_columns)['seg_genbank'].apply(lambda v: len(v))

# merge them back to one DF
new_refseq =  pandas.merge(x_new_refseq, x_new_refseq_ct,  on=groupby_columns, how='outer', indicator=True)
new_refseq.rename(columns = {'seg_refseq_x':'seg_refseq', 'seg_refseq_y':'seq_refseq_ct', '_merge':'seg_refseq_merge'}, inplace= True)

new_genbank = pandas.merge(x_new_genbank,x_new_genbank_ct, on=groupby_columns, how='outer', indicator=True)
new_genbank.rename(columns = {'seg_genbank_x':'seg_genbank', 'seg_genbank_y':'seq_genbank_ct', '_merge':'seg_genbank_merge'}, inplace= True)

tabulated_df = pandas.merge(new_genbank,new_refseq,        on=groupby_columns, how='outer', indicator=True)

# build SQL
def generate_update_sql(r):
    ea=r['Exemplar_Additional'].upper()
    if ea=='E' or ea=='A':
        # common columns
        ea_col             ='exemplar'
        refseq_org_col     ='refseq_organism' # not implemented, yet
        refseq_taxids_col  ='refseq_taxids' # not implemented, yet
        # SQL prefix
        sql_prefix= ""
        # select columns by type
        if   ea =='E':
            genbank_acc_col ='exemplar_genbank_accession'
            refseq_acc_col  ='exemplar_refseq_accession'
        elif ea =='A':
            genbank_acc_col ='isolate_genbank_accession_csv'
            refseq_acc_col  ='isolate_refseq_accession'  # just added
    else:
        # error message as SQL comment
        sql_prefix = "-- ERROR: Exemplar_Additional='"+ea+"' for "

    # format SQL statement
    sql = (
        sql_prefix +
        "UPDATE [vmr] SET " +
        "  " + genbank_acc_col + "= '" + r['seg_genbank'] +"' "
        ", " + refseq_acc_col  + "= '" + r['seg_refseq']  +"' "
        " FROM [vmr] " +
        " WHERE [" + genbank_acc_col + "] = '" + r['Original_Accession_String']+"' "
        " AND   [" + ea_col           + "] = '" + ea +"' "
        " AND   [" +'species'        + "] = '" + r['Species'] +"' "
    )
    # final result
    return(sql)

# generate sql update statements for all rows
tabulated_df['sql_update'] = tabulated_df.apply(lambda r: generate_update_sql(r), axis=1)

tabulated_fname="out.tsv"
if args.verbose: print("Writing: ", tabulated_fname)
pandas.DataFrame.to_csv(tabulated_df, tabulated_fname, sep="\t")
if args.verbose: print("   Write {0} rows, {1} columns.".format(*tabulated_df.shape))
