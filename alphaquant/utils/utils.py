import os

def get_condpairname(condpair):
    return f"{condpair[0]}_VS_{condpair[1]}"

def get_condpair_from_condpairname(condpairname):
    return condpairname.split("_VS_")


from anytree.importer import JsonImporter
def read_condpair_tree(cond1, cond2, results_folder = os.path.join(".", "results")):
    """reads the merged and clustered iontree for a given condpair"""
    condpairname = get_condpairname([cond1, cond2])
    tree_file =os.path.join(results_folder, f"{condpairname}.iontrees.json")
    if not os.path.isfile(tree_file):
        return None
    importer = JsonImporter()
    filehandle = open(tree_file, 'r')
    jsontree = importer.read(filehandle)
    filehandle.close()
    return jsontree


def cut_trailing_parts_seqstring(seqstring):
    return seqstring.replace("SEQ_", "").rstrip("_")