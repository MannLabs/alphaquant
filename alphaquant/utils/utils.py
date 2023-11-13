import os
from anytree.importer import JsonImporter
import glob


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



def cut_trailing_parts_seqstring(seqstring):
    return seqstring.replace("SEQ_", "").rstrip("_")


def get_condpairname(condpair):
    return f"{condpair[0]}_VS_{condpair[1]}"

def get_condpair_from_condpairname(condpairname):
    return condpairname.split("_VS_")