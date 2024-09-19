#!/usr/bin/env python3
#
# VMR accession list parsing
#
print("accession_utils.py:import...")
import sys
import pandas as pd
from natsort import natsorted

au_tmi = False

def parse_seg_accession_list(isolate_id,acc_list_str):
    # don't loose it on nan's
    if pd.isna(acc_list_str):
        return {}
        #acc_list_str = ""

    # remove whitespace.
    acc_list_str = acc_list_str.replace(" ","")

    # instead of trying to split by commas and semicolons, I just replace the commas with semicolons. 
    acc_list_str = acc_list_str.replace(",",";")

    # split into list: ";" 
    accession_list = acc_list_str.split(';')
    if au_tmi: print("accession_list:"+"|".join(accession_list))

    # 
    # for each [SEG:]ACCESSION
    # 
    sa_dict = {} # seg_name-accession map
    accession_index = 0
    for seg_acc_str in accession_list:
        if au_tmi: print("seg_acc_str:"+seg_acc_str)

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
                if au_tmi: print("sa_dict["+str(accession_index)+"]:"+str(sa_dict[accession_index]))
            elif len(seg_acc_pair)==2:
                # seg_name:accession
                segment_name = seg_acc_pair[0]
                accession = seg_acc_pair[1]
                sa_dict[segment_name] = {"accession":accession, "segment_name":segment_name, "accession_index":accession_index, "isolate_id":isolate_id}
                if au_tmi: print("sa_dict["+str(segment_name)+"]:"+str(sa_dict[segment_name]))

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
    if len(gb_accessions_dict) != len(rs_accessions_dict) and len(rs_accessions_dict)>0:
        print("WARNING[isolate:"+str(isolate_id)+"]: gb_n_acc: "+str(len(gb_accessions_dict))+" != rs_n_acc:"+str(len(rs_accessions_dict)),file=sys.stderr)
        #print("gb_accessions_dict:\n",gb_accessions_dict)
        #print("rs_accessions_dict:\n",rs_accessions_dict)
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
                'segment_name': gb_entry.get('segment_name', rs_entry.get('segment_name')),
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
                'segment_name': gb_entry.get('segment_name', rs_entry.get('segment_name')),
                'gb_segment_name': gb_entry.get('segment_name', None),
                'gb_accession'   : gb_entry.get('accession', None),
                'rs_segment_name': rs_entry.get('segment_name', None),
                'rs_accession'   : rs_entry.get('accession', None)
            }                

    #convert back to an ordered list/data frame
    #print("------------\n ## merged_dict:\n------------\n", merged_dict)
    values = list(merged_dict.values())
    sorted_values = natsorted(values, key=lambda x: [x['accession_index'],x['segment_name']])

    # Print the entries in the sorted order
    merged_df = pd.DataFrame(sorted_values)
    #print("------------\n ## merged_df:\n------------\n", merged_df)

    return merged_df

