import os
from anytree.importer import JsonImporter
import glob
import pandas as pd
import re

NODETYPE2REGEX ={'frgion': 'SEQ.*MOD.*CHARGE.*FRG)(ION.*)',
 'ms1_isotopes': '(SEQ.*MOD.*CHARGE.*MS1)(ISO.*)',
 'mod_seq_charge': 'SEQ.*MOD.*CHARGE.*)(FRG.*|MS1.*)',
 'mod_seq': 'SEQ.*MOD.*)(CHARGE.*)',
 'seq': 'SEQ.*)(MOD.*)'}


def read_all_trees_in_results_folder(results_folder):
    condpair2tree = {}

    for json_file in glob.glob(os.path.join(results_folder, '*.iontrees.json')):
        condpair_tree = read_tree_from_json(json_file)
        condpairname = get_condpairname(condpair_tree.name) #name is a tuple of cond1 and cond2
        condpair2tree[condpairname] = condpair_tree

    return condpair2tree

def read_condpair_tree(cond1, cond2, results_folder = os.path.join(".", "results")):
    """reads the merged and clustered iontree for a given condpair"""
    condpairname = get_condpairname([cond1, cond2])
    tree_file =os.path.join(results_folder, f"{condpairname}.iontrees.json")
    if not os.path.isfile(tree_file):
        return None
    
    return read_tree_from_json(tree_file)

def read_tree_from_json(tree_file):
    importer = JsonImporter()
    filehandle = open(tree_file, 'r')
    jsontree = importer.read(filehandle)
    filehandle.close()
    return jsontree

def read_all_results_files_in_results_folder(results_folder):
    condpair2df = {}

    for results_file in glob.glob(os.path.join(results_folder, '*.results.tsv')):
        results_df = pd.read_csv(results_file, sep = "\t")
        results_filename_nopath = os.path.basename(results_file)
        condpairname = results_filename_nopath.replace(".results.tsv", "")
        condpair2df[condpairname] = results_df

    return condpair2df

def cut_trailing_parts_seqstring(seqstring):
    return seqstring.replace("SEQ_", "").rstrip("_")

def get_condpairname(condpair):
    return f"{condpair[0]}_VS_{condpair[1]}"

def get_condpair_from_condpairname(condpairname):
    return condpairname.split("_VS_")

def remove_file_extension(filename):
    trimmed_filename = os.path.splitext(filename)[0]
    return trimmed_filename

def convert_ion_string_to_node_type(ionstring, node_type): #for example I have a full quant_id that describes a fragment ion, I want to shorten it to the specified leve, e.g. sequence
    regex = NODETYPE2REGEX[node_type]
    match = re.match(regex, ionstring)
    if match:
        return match.group(1)
    else:
        raise ValueError(f"Could not match {ionstring} to {node_type}. This function only works for the following node types: seq, mod_seq, mod_seq_charge")


