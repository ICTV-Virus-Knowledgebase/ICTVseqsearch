#!/usr/bin/env python3
#
# test VMR accession list parsing
#
print("import...")
import sys
import pandas as pd
from natsort import natsorted

# our local package for "[Segment:]Accession[;]" list parsing
sys.path.append('..')
from accession_utils import parse_seg_accession_list,merge_acc_dicts

def print_dict_order_by_accession_index(merged_dict):
    # Sort the dictionary items by 'accession_index' key
    sorted_items = sorted(merged_dict.items(), key=lambda item: item[1]['accession_index'])

    # Print the entries in the sorted order
    for key, value in sorted_items:
        print(f"Key: {key}, Entry: {value}")


def write_dict_to_tsv(merged_dict, file_name="parse_seg_acc_list.accessions_map.tsv"):
    # Convert the merged dictionary into a pandas DataFrame
    df = pd.DataFrame.from_dict(merged_dict, orient='index')

    # Sort the DataFrame by 'isolate_id' and 'accession_index'
    df = df.sort_values(by=['isolate_id', 'accession_index','segment_name'])

    # Write the DataFrame to a TSV file
    df.to_csv(file_name, sep='\t', index=False) #index_label='key')

data = [
    { "case":"acc_both_empty", "isolate_id":1003732.1, "gb_a":'', "rs_a":'' },
    { "case":"acc_rs_missing", "isolate_id":1003732.2, "gb_a":'HM246720', "rs_a":'' },
    { "case":"acc_gb_missing", "isolate_id":1003732.3, "gb_a":'', "rs_a":'NC_027989' },
    { "case":"acc_single",     "isolate_id":1003732.4, "gb_a":'HM246720', "rs_a":'NC_027989' },
    { "case":"acc_list",       "isolate_id":1003732.5, "gb_a":'HM246720; HM246721; HM246722; HM246723; HM246724', "rs_a":'NC_027989; NC_041833; NC_041831; NC_041832; NC_041834' },
    { "case":"seg_acc_single_matched",    "isolate_id":1007556.1, "gb_a":'DNA-C: EF546812', "rs_a":'DNA-C: NC_010318'},
    { "case":"seg_acc_single_mismatched", "isolate_id":1007556.2, "gb_a":'DNA-C: EF546812', "rs_a":'DNA-N: NC_010318'},
    { "case":"seg_acc_single_ragged",     "isolate_id":1007556.3, "gb_a":'DNA-C: EF546812', "rs_a":'DNA-C: NC_010318 ;  DNA-U3:     NC_010315'},
    { "case":"seg_acc_list",              "isolate_id":1007556.4, "gb_a":'DNA-C: EF546812; DNA-M: EF546811; DNA-N: EF546808; DNA-R: EF546813; DNA-S:EF546810; DNA-U3: EF546809', "rs_a":'DNA-C: NC_010318; DNA-M: NC_010317; DNA-N:NC_010314; DNA-R: NC_010319; DNA-S: NC_010316; DNA-U3:     NC_010315'}
]

meta_df = pd.DataFrame()

for datum in data:
    print("# ----------------------------------------------------------------------")
    print("# test case: "+datum["case"])
    print("# ----------------------------------------------------------------------")

    gb_d = parse_seg_accession_list(datum["isolate_id"],datum["gb_a"])
    print("# --- GENBANK: ---\n\t"); print_dict_order_by_accession_index(gb_d) # ,str(gb_d))
    rs_d = parse_seg_accession_list(datum["isolate_id"],datum["rs_a"])
    print("# --- REFSEQ:  ---\n\t"); print_dict_order_by_accession_index(rs_d) #,str(rs_d))

    merged_df = merge_acc_dicts(datum["isolate_id"],gb_d,rs_d)
    print("# --- MERGE:   ---\n")
    fname="parse_seg_acc_list."+datum["case"]+".tsv"
    merged_df.to_csv(fname, sep='\t', index=False)
#    #print_dict_order_by_accession_index(merge_d)
#    write_dict_to_tsv(merge_d, fname)
    print("Wrote "+fname)

    merged_df['case'] = datum["case"]
    #print("------------\n ## merged_df: w/ case \n------------\n", merged_df)
    
    meta_df = pd.concat([meta_df, merged_df], ignore_index=True)
    #print("------------\n ## meta_df: w/ case \n------------\n", meta_df)
    #meta_df.extend(merged_df) # extend is for dict's

# final results to DF
#meta_df = pd.DataFrame(meta_list)

# Sort the DataFrame by 'isolate_id' and 'accession_index'
print("meta_df.columns:",str(meta_df.columns))
print(meta_df)
meta_df= meta_df.sort_values(by=['case','isolate_id', 'accession_index','segment_name'])


# Write the DataFrame to a TSV file
file_name = "parse_seg_acc_list.meta_list.tsv"
meta_df.to_csv(file_name, sep='\t', index=False) #, index_label='key')
