#!/usr/bin/env python3
#
# Generate VMR RefSeq strings, and SQL updates to load them in the db
#
# INPUT: vmr.tsv
#    ID', 'Sort', 'Isolate Sort', 'Realm', ..., 'Species', 'Exemplar or additional isolate', 'Virus name(s)', 'Virus name abbreviation(s)',
#    'Virus isolate designation', 'Virus GENBANK accession', 'Virus REFSEQ accession', 'Genome coverage', 'Genome composition', 'Host source'
# INPUT: processed_accessions_b.tsv
#    Species Accession_IDs segment Sort Original_Accession_String
# INPUT: genbank_refseq_map_wide_fixed.txt
#    REFSEQ GENBANK NCBI_ISOLATE_NAME ...   
#
# Update will use Original_Accession_String as the PK into VMR table.
# It will update both the Genbank and RefSeq accession strings
#
print("Importing...")
import argparse
import pandas
from natsort import natsort_keygen
from pprint import pprint
import re
print("Parsing...")
parser = argparse.ArgumentParser(description="")

#setting arguments.
parser.add_argument('-verbose',help="printout details during run",action=argparse.BooleanOptionalAction, default=True)
parser.add_argument("-vmr",help="TSV: VMR xlsx translated to .tsv", default="vmr.tsv")
parser.add_argument("-vmr_accessions",help="TSV: VMR accessions, one per line, with segment name", default="processed_accessions_b.tsv")
parser.add_argument("-refseq_map",help="TSV: maps Genbank accessios to matching RefSeq accessions", default="genbank_refseq_map_wide_fixed2.txt")
parser.add_argument("-sql",help="output filename for update SQL.",default="vmr_genbank_refseq_update.sql")
parser.add_argument("-tsv",help="output filename for updated refseq TSV", default="vmr_genbank_refseq_update.tsv")
args = parser.parse_args()

print("###############################################################################################################")
print("# Loads excel inputs")
print("###############################################################################################################")

if args.verbose: print("Loading: ", args.vmr)
vmr_df = pandas.read_csv(args.vmr, sep="\t")
if args.verbose: print("   Read {0} rows, {1} columns.".format(*vmr_df.shape))
if args.verbose: print("   columns:", list(vmr_df.columns))
vmr_df.rename(columns = {'Isolate ID':'Isolate_ID'}, inplace= True)
vmr_df['Virus NCBI taxid'] = '' #list(map(lambda v: '',vmr_acc_df['Original_GENBANK_Accessions']))

if args.verbose: print("Loading: ", args.vmr_accessions)
vmr_acc_df = pandas.read_csv(args.vmr_accessions, sep="\t")
if args.verbose: print("   Read {0} rows, {1} columns.".format(*vmr_acc_df.shape))
if args.verbose: print("   columns:", list(vmr_acc_df.columns))
# fake an original field for ncbi_taxid
vmr_acc_df['Virus NCBI taxid'] = '' #list(map(lambda v: '',vmr_acc_df['Original_GENBANK_Accessions']))
if args.verbose: print("   columns:", list(vmr_acc_df.columns))

if args.verbose: print("Loading: ", args.refseq_map)
# name "refseq" Accession_IDs, so column names match for merge
refseq_colnames=['refseq', 'Accession_IDs', 'ncbi_virus_name']
# genbank_refseq_map_wide.txt
print("check: args.refseq_map="+args.refseq_map+" for wide gives",re.match(r".*wide.*",args.refseq_map))
if bool(re.match(r".*wide.*",args.refseq_map)):
    refseq_colnames=['refseq', 'genbank', 'ncbi_virus_name','ncbi_taxid','ncbi_slen','ncbi_moltype','ncbi_topology','ncbi_completeness','ncbi_strand']

refseq_map = pandas.read_csv(args.refseq_map, dtype=str, sep="\t", header=None, names=refseq_colnames)
if args.verbose: print("   Read {0} rows, {1} columns.".format(*refseq_map.shape))
if args.verbose: print("   columns:", list(refseq_map.columns))
if args.verbose: refseq_map.info()

print("###############################################################################################################")
print("# clean data ")
print("###############################################################################################################")

#
# remove .# suffix from accession numbers
#
vmr_acc_df[      "genbank"] = vmr_acc_df.Accession.replace(    to_replace='\\.[-0-9]+$',  regex=True, value='')
vmr_acc_df[      "genbank"] = vmr_acc_df.genbank.replace(      to_replace='\\([.0-9]*\\)', regex=True, value='')
refseq_map[  "genbank"] = refseq_map.genbank.replace(  to_replace='\\.[-0-9]+$',  regex=True, value='')

print("###############################################################################################################")
print("# merge data frames")
print("###############################################################################################################")

vmr_refseq_merge=pandas.merge(vmr_acc_df, refseq_map, on=['genbank'], how="outer", indicator=True)

vmr_refseq_raw=vmr_refseq_merge.query('_merge!="right_only"')
vmr_refseq_bad=vmr_refseq_merge.query('_merge=="right_only"')

bad_fname = "vmr_genbank_refseq_unmatched.tsv"
if args.verbose: print("Writing: ", bad_fname)
pandas.DataFrame.to_csv(vmr_refseq_bad, bad_fname, sep="\t")
if args.verbose: print("   Wrote {0} rows, {1} columns.".format(*vmr_refseq_bad.shape))

if args.verbose: print("Writing: ", args.tsv)
pandas.DataFrame.to_csv(vmr_refseq_raw, args.tsv, sep="\t")
if args.verbose: print("   Wrote {0} rows, {1} columns.".format(*vmr_refseq_raw.shape))
if args.verbose: print("   Columns:", list(vmr_refseq_raw.columns))

# filter to only records with refseq (hum - what about removing refseqs?)
#if args.verbose: print("Filter out VMR entries without RefSeq matches")
#vmr_refseq_map = vmr_refseq_raw.dropna(subset=['refseq']).copy()
#pandas.DataFrame.to_csv(vmr_refseq_map, args.tsv, sep="\t")
#if args.verbose: print("   Wrote {0} rows, {1} columns.".format(*vmr_refseq_map.shape))
vmr_refseq_map = vmr_refseq_raw

print("###############################################################################################################")
print("# label accessions/taxids with segment names")
print("###############################################################################################################")


# Define a function to concatenate segment and refseq values
def format_seg_acc(seg,acc):
    if seg != seg: 
        if acc != acc:
            return ''
        else:
            return str(acc)
    if acc != acc:
        return str(seg)+': '
    else:
        return str(seg)+': '+str(acc)

# Format accessions with segment names
vmr_refseq_map['seg_genbank']=    vmr_refseq_map.apply(lambda r: format_seg_acc(r['Segment_Name'],r['genbank']),axis=1)
vmr_refseq_map['seg_refseq']=     vmr_refseq_map.apply(lambda r: format_seg_acc(r['Segment_Name'],r['refseq']), axis=1)
vmr_refseq_map['seg_ncbi_taxid']= vmr_refseq_map.apply(lambda r: format_seg_acc(r['Segment_Name'],r['ncbi_taxid']), axis=1)


print("###############################################################################################################")
print("# tabulate segments into lists")
print("###############################################################################################################")


##
##
## Group by Isolate_ID, Sort by Accession Name/Idx and apply the concatenation function
##
##
groupby_columns = ['Isolate_ID']

##
print("# GENBANK")
##

print("# count genbank segments")
new_genbank_ct= vmr_refseq_map[[*groupby_columns,'seg_genbank']].groupby(by=groupby_columns)['seg_genbank'].apply(lambda v: len(v)).reset_index()
# rename output column
new_genbank_ct.rename(columns = {'seg_genbank':'seg_genbank_ct'}, inplace= True)
print("new_genbank_ct: {0} rows, {1} columns.".format(*new_genbank_ct.shape))
# merge into VMR
vmr_out_df = pandas.merge(vmr_df,new_genbank_ct, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))
#print("vmr_out_df + :",list(vmr_out_df.columns))

print("# merge seg:genbank into list (orig order)")
new_genbank_list= vmr_refseq_map[[*groupby_columns,'seg_genbank','Accession_Index']].groupby(by=groupby_columns)[['seg_genbank','Accession_Index']].apply(lambda v: '; '.join(v.sort_values(by="Accession_Index")['seg_genbank'])).reset_index()
# rename output column
new_genbank_list.rename(columns = {0:'seg_genbank_idx_list'}, inplace= True)
print("new_genbank_list: {0} rows, {1} columns.".format(*new_genbank_list.shape))
#print("new_genbank_list:",list(new_genbank.columns))
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,new_genbank_list, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))

print("# merge seg:genbank into list (seg name natural order, then acc_index))")
nat_genbank_list= vmr_refseq_map[[*groupby_columns,'seg_genbank','Segment_Name','Accession_Index']].groupby(by=groupby_columns)[['seg_genbank','Accession_Index','Segment_Name']].apply(lambda v: '; '.join(v.sort_values(by=["Segment_Name","Accession_Index"],key=natsort_keygen())["seg_genbank"])).reset_index()
# rename output column
nat_genbank_list.rename(columns = {0:'seg_genbank_nat_list'}, inplace= True)
print("nat_genbank_list: {0} rows, {1} columns.".format(*nat_genbank_list.shape))
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,nat_genbank_list, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))

##
print("# REFSEQ")
##      

print("# count uniq refseq accessions")
new_refseq_ct_uniq= vmr_refseq_map[[*groupby_columns,'refseq']].groupby(by=groupby_columns)['refseq'].apply(lambda v: len(set(v)-{''})).reset_index()
# rename output column
new_refseq_ct_uniq.rename(columns = {'refseq':'refseq_ct_uniq'}, inplace= True)
print("new_refseq_ct_uniq: {0} rows, {1} columns.".format(*new_refseq_ct_uniq.shape))
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,new_refseq_ct_uniq, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))
print("vmr_out_df + :",list(vmr_out_df.columns))

print("# first refseq accession")
new_refseq_first= vmr_refseq_map[[*groupby_columns,'refseq']].groupby(by=groupby_columns)['refseq'].apply(lambda v: next(iter(v),"")).reset_index()
# rename output column
new_refseq_first.rename(columns = {'refseq':'refseq_first'}, inplace= True)
print("new_refseq_first: {0} rows, {1} columns.".format(*new_refseq_first.shape))
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,new_refseq_first, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))
#print("vmr_out_df + :",list(vmr_out_df.columns))

print("# merge seg:refseq into list (orig order)")
new_refseq=vmr_refseq_map[[*groupby_columns,'seg_refseq','Accession_Index']].groupby(by=groupby_columns)[['seg_refseq','Accession_Index']].apply(lambda v: '; '.join(v.sort_values(by="Accession_Index")['seg_refseq'])).reset_index()
new_refseq.rename(columns = {0:'seg_refseq_idx_list'}, inplace= True)
print("new_refseq:",list(new_refseq.columns))
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,new_refseq, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))

print("# merge seg:genbank into list (seg name natural order, then acc_index))")
nat_refseq=vmr_refseq_map[[*groupby_columns,'seg_refseq','Segment_Name','Accession_Index']].groupby(by=groupby_columns)[['seg_refseq','Accession_Index','Segment_Name' ]].apply(lambda v: '; '.join(v.sort_values(by=["Segment_Name","Accession_Index"],key=natsort_keygen())["seg_refseq"])).reset_index()
nat_refseq.rename(columns = {0:'seg_refseq_nat_list'}, inplace= True)
print("nat_refseq:",list(nat_refseq.columns))
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,nat_refseq, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))

# chose '' if all accessions are '', else choose the list
def select_refseq_empty_or_list(uct,acc,seg_acc_list):
    if uct > 1:
        return seg_acc_list
    if  acc=='' or acc!=acc:
        return ''
    return seg_acc_list

print("# create seq_refseq_nat/idx_empty")
vmr_out_df['seg_refseq_idx_list_empty']= vmr_out_df.apply(lambda r: select_refseq_empty_or_list(r['refseq_ct_uniq'],r['refseq_first'],r['seg_refseq_idx_list']), axis=1)
vmr_out_df['seg_refseq_nat_list_empty']= vmr_out_df.apply(lambda r: select_refseq_empty_or_list(r['refseq_ct_uniq'],r['refseq_first'],r['seg_refseq_nat_list']), axis=1)
if args.verbose: print("   columns:", list(vmr_out_df.columns))

##
##
print("# TAXID") # ----------------------------------------------------------------------
##
## 

print("# count number of uniq ncbi_taxid per isolate (some segments disagree!)")
taxid_ct_uniq= vmr_refseq_map[[*groupby_columns,'ncbi_taxid']].groupby(by=groupby_columns)['ncbi_taxid'].apply(lambda v:  len(set(v)-{''})).reset_index()
# rename output column
taxid_ct_uniq.rename(columns = {'ncbi_taxid':'ncbi_taxid_ct_uniq'}, inplace= True)
print("taxid_ct_uniq: {0} rows, ? columns.".format(*taxid_ct_uniq.shape))
print(taxid_ct_uniq[0:1])
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,taxid_ct_uniq, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))

print("# first ncbi_taxid")
taxid_first= vmr_refseq_map[[*groupby_columns,'ncbi_taxid']].groupby(by=groupby_columns)['ncbi_taxid'].apply(lambda v: next(iter(v), "")).reset_index()
# rename output column
taxid_first.rename(columns = {'ncbi_taxid':'ncbi_taxid_first'}, inplace= True)
print("taxid_first: {0} rows, ? columns.".format(*taxid_first.shape))
print(taxid_first[0:1])
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,taxid_first, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))

print("# merge seg_taxid into list by accession_index")
taxid_idx_list=vmr_refseq_map[[*groupby_columns,'seg_ncbi_taxid','Accession_Index']].groupby(by=groupby_columns)[['seg_ncbi_taxid','Accession_Index']].apply(lambda v: '; '.join([str(element) for element in v.sort_values(by="Accession_Index")['seg_ncbi_taxid']])).reset_index()
# rename output column
taxid_idx_list.rename(columns = {0:'seg_ncbi_taxid_idx_list'}, inplace= True)
print("taxid_idx_list:",list(taxid_idx_list.columns))
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,taxid_idx_list, on=groupby_columns, how='outer', indicator=False)

print("# merge seg_taxid into list by (natural segment name, accession_index)")
taxid_nat_list=vmr_refseq_map[[*groupby_columns,'seg_ncbi_taxid','Segment_Name','Accession_Index']].groupby(by=groupby_columns)[['seg_ncbi_taxid','Accession_Index','Segment_Name' ]].apply(lambda v: '; '.join([str(element) for element in v.sort_values(by=["Segment_Name","Accession_Index"],key=natsort_keygen())["seg_ncbi_taxid"]])).reset_index()
# rename output column
taxid_nat_list.rename(columns = {0:'seg_ncbi_taxid_nat_list'}, inplace= True)
print("taxid_nat_list:",list(taxid_nat_list.columns))
# merge into VMR
vmr_out_df = pandas.merge(vmr_out_df,taxid_nat_list, on=groupby_columns, how='outer', indicator=False)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))

# chose singluar taxid (1 seg, or all segments agree) or "seg:taxid;" list, if segments
# discordant taxids
# Define a function to concatenate segment and refseq values
def select_taxid_uniq_or_list(uct,taxid,taxid_list):
    if uct > 1:
        return taxid_list
    return taxid
print("# create seq_ncbi_taxid_new/nat_uniq")
if args.verbose: print("   columns:", list(vmr_out_df.columns))

vmr_out_df['seg_ncbi_taxid_idx_uniq']= vmr_out_df.apply(lambda r: select_taxid_uniq_or_list(r['ncbi_taxid_ct_uniq'],r['ncbi_taxid_first'],r['seg_ncbi_taxid_idx_list']), axis=1)
vmr_out_df['seg_ncbi_taxid_nat_uniq']= vmr_out_df.apply(lambda r: select_taxid_uniq_or_list(r['ncbi_taxid_ct_uniq'],r['ncbi_taxid_first'],r['seg_ncbi_taxid_nat_list']), axis=1)
print("vmr_out_df: {0} rows, {1} columns.".format(*vmr_out_df.shape))
print("")


print("#")
print("# final vmr_out_df ")
print("#")
if args.verbose: print("   Merged: {0} rows, {1} columns.".format(*vmr_out_df.shape))
if args.verbose: print("   columns:", list(vmr_out_df.columns))
print(vmr_out_df[0:1])
print(" ")


# build SQL
#  columns                field (index order)       field (natural seg, idx order)   orig_field
#  genbank_accessions     seg_genbank_idx_list      seg_genbank_nat_list             'Virus GENBANK accession'
#  refseq_accessions      seg_refseq_idx_list       seg_refseq_nat_list              'Virus REFSEQ accesssion'
#  refseq_taxids          seg_ncbi_taxid_idx_uniq   seg_ncbi_taxid_nat_uniq          'Virus NCBI taxid'
def generate_update_sql(r,col,orig_field,field):
    #print("ROW: ",r)

    # format SQL statement
    # sue the CASE to make '' into NULL in SQL
    sql = (
        "UPDATE [species_isolates] SET " +
        "  ["+col+"] = (case when '"+str(r[field])+"'<>'' then '"+str(r[field])+"' end) "
        " FROM [species_isolates] " +
        " WHERE [isolate_id] = " + str(int(r['Isolate_ID'])) + " " 
        " AND ( "
        " ( ["+col+"] <> '" + str(r[field]) +"' "+ " AND  ["+col+"] NOT LIKE '%(%.%)') "   
        " OR (["+col+"] is NULL AND '" + str(r[field]) +"'<>'') "
        ")"
        " -- ["+orig_field+"]='"+str(r[orig_field])+"'"
    )
    # final result
    return(sql if r[orig_field]!=r[field] else '')

# generate sql update statements for all rows
vmr_out_df['sql_update_genbank'] = vmr_out_df.apply(lambda r: generate_update_sql(r,'genbank_accessions','Virus GENBANK accession',  'seg_genbank_nat_list'       ), axis=1)
vmr_out_df['sql_update_refseq']  = vmr_out_df.apply(lambda r: generate_update_sql(r,'refseq_accessions', 'Virus REFSEQ accession',   'seg_refseq_nat_list_empty'  ), axis=1)
vmr_out_df['sql_update_taxid']   = vmr_out_df.apply(lambda r: generate_update_sql(r,'refseq_taxids',     'Virus NCBI taxid',         'seg_ncbi_taxid_nat_uniq'    ), axis=1)

vmr_out_fname="out.tsv"
if args.verbose: print("Writing: ", vmr_out_fname)
pandas.DataFrame.to_csv(vmr_out_df, vmr_out_fname, sep="\t")
if args.verbose: print("   Write {0} rows, {1} columns.".format(*vmr_out_df.shape))
