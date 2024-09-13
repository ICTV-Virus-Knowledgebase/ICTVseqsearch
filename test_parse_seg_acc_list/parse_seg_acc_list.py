#!/usr/bin/env python3
#
# test VMR accession list parsing
#
print("import...")
import sys
import pandas as pd

def parse_seg_accession_list(isolate_id,acc_list_str):
    # remove whitespace.
    acc_list_str = acc_list_str.replace(" ","")

    # instead of trying to split by commas and semicolons, I just replace the commas with semicolons. 
    acc_list_str = acc_list_str.replace(",",";")

    # split into list: ";" 
    accession_list = acc_list_str.split(';')
    #if True or args.tmi: print("accession_list:"+"|".join(accession_list))

    # 
    # for each [SEG:]ACCESSION
    # 
    sa_dict = {} # seg_name-accession map
    accession_index = 0
    for seg_acc_str in accession_list:
        #if True or args.tmi: print("seg_acc_str:"+seg_acc_str)

        # track accession/segment order, so it can be preserved
        accession_index += 1

        # split optional "segment_name:" prefix on accessions
        seg_acc_pair = seg_acc_str.split(':')
        segment_name = None
        accession    = None
        if len(seg_acc_pair)==0 or len(seg_acc_pair)>2:
            print("ERROR[isolate_id:"+str(isolate_id)+": [seg:]acc >1 colon: '"+str(seg_acc_pair)+"' from '"+acc_list_str+"'",file=sys.stderr)
        else:
            if len(seg_acc_pair)==1:
                # bare accession
                accession = seg_acc_pair[0]
                sa_dict[accession_index] = {"accession":accession, "segment_name":None, "accession_index":accession_index, "isolate_id":isolate_id}
                #if True or args.tmi: print("sa_dict["+str(accession_index)+"]:"+str(sa_dict[accession_index]))
            elif len(seg_acc_pair)==2:
                # seg_name:accession
                segment_name = seg_acc_pair[0]
                accession = seg_acc_pair[1]
                sa_dict[segment_name] = {"accession":accession, "segment_name":segment_name, "accession_index":accession_index, "isolate_id":isolate_id}
                #if True or args.tmi: print("sa_dict["+str(segment_name)+"]:"+str(sa_dict[segment_name]))

            # QC accessions
            number_count = 0
            letter_count = 0
            # counting letters
            for char in accession:
                if char in 'qwertyuiopasdfghjklzxcvbnm':
                    letter_count = letter_count+1
            # counting numbers
                elif char in '1234567890':
                    number_count = number_count+1
            #checks if current selection fits what an accession number should be
            if accession != "" and not ((len(str(accession)) == 8 or 6) and letter_count<3 and number_count>3):
                print("ERROR[isolate_id:"+str(isolate_id)+"]: suspect accesssion '"+accession+"'",file=sys.stderr)

                
    # we'll check later if this segment has a name 
    return(sa_dict)

def print_dict_order_by_accession_index(merged_dict):
    # Sort the dictionary items by 'accession_index' key
    sorted_items = sorted(merged_dict.items(), key=lambda item: item[1]['accession_index'])

    # Print the entries in the sorted order
    for key, value in sorted_items:
        print(f"Key: {key}, Entry: {value}")

def merge_acc_dicts(isolate_id, gb_accessions_dict, rs_accessions_dict):
    # merge parallel lists (not nice)
    if len(gb_accessions_dict) != len(rs_accessions_dict):
        print("WARNING[isolate:"+str(isolate_id)+"]: gb_n_acc: "+str(len(gb_accessions_dict))+" != rs_n_acc:"+str(len(rs_accessions_dict)),file=sys.stderr)

    # iterate over gb, build merged list
    merged_dict={}

    # Get the union of keys from both dictionaries
    all_keys = set(gb_accessions_dict.keys()).union(rs_accessions_dict.keys())

    # Loop over each key in the union
    for key in all_keys:
        # Get the corresponding dictionary for the key from each dictionary
        gb_entry = gb_accessions_dict.get(key, {})
        rs_entry = rs_accessions_dict.get(key, {})

        # Merge the two dictionaries for the given key
        if issubclass(type(key),int):
            # index (no segment names)
            merged_dict[key] = {
                'isolate_id': gb_entry.get('isolate_id', rs_entry.get('isolate_id')),
                'accession_index': gb_entry.get('accession_index', rs_entry.get('accession_index')),
                'gb_segment_name': None,
                'gb_accession': gb_entry.get('accession', None),
                'rs_segment_name': None,
                'rs_accession': rs_entry.get('accession', None)
            }                
        else:
            # seg name
            merged_dict[key] = {
                'isolate_id': gb_entry.get('isolate_id', rs_entry.get('isolate_id')),
                'accession_index': gb_entry.get('accession_index', rs_entry.get('accession_index')),
                'gb_segment_name': gb_entry.get('segment_name', None),
                'gb_accession'   : gb_entry.get('accession', None),
                'rs_segment_name': rs_entry.get('segment_name', None),
                'rs_accession'   : rs_entry.get('accession', None)
            }                

    return merged_dict        # merge parallel lists (not nice)

def write_dict_to_tsv(merged_dict, file_name="parse_seg_acc_list.accessions_map.tsv"):
    # Convert the merged dictionary into a pandas DataFrame
    df = pd.DataFrame.from_dict(merged_dict, orient='index')

    # Sort the DataFrame by 'isolate_id' and 'accession_index'
    df = df.sort_values(by=['isolate_id', 'accession_index'])

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

meta_list = []

for datum in data:
    print("# ----------------------------------------------------------------------")
    print("# test case: "+datum["case"])
    print("# ----------------------------------------------------------------------")

    gb_d = parse_seg_accession_list(datum["isolate_id"],datum["gb_a"])
    print("# --- GENBANK: ---\n\t"); print_dict_order_by_accession_index(gb_d) # ,str(gb_d))
    rs_d = parse_seg_accession_list(datum["isolate_id"],datum["rs_a"])
    print("# --- REFSEQ:  ---\n\t"); print_dict_order_by_accession_index(rs_d) #,str(rs_d))

    merge_d = merge_acc_dicts(datum["isolate_id"],gb_d,rs_d)
    print("# --- MERGE:   ---\n")
    print_dict_order_by_accession_index(merge_d)
    fname="parse_seg_acc_list."+datum["case"]+".tsv"
    write_dict_to_tsv(merge_d, fname)
    print("Wrote "+fname)

    for item in merge_d.values(): 
        item.update({'case': datum["case"]})

    meta_list.extend(merge_d.values())

# final results to DF
meta_df = pd.DataFrame(meta_list)

# Sort the DataFrame by 'isolate_id' and 'accession_index'
meta_df= meta_df.sort_values(by=['case','isolate_id', 'accession_index'])


# Write the DataFrame to a TSV file
file_name = "parse_seg_acc_list.meta_list.tsv"
meta_df.to_csv(file_name, sep='\t', index=False) #, index_label='key')
