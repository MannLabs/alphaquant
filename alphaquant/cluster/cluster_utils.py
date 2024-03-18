
# Cell
import anytree
from statistics import NormalDist
import numpy as np
import collections
import alphaquant.config.variables as aqvariables
from anytree import Node, LevelOrderGroupIter

TYPES = ["base","frgion", "ms1_isotopes", "mod_seq_charge", "mod_seq", "seq", "gene"]
LEVELS = ["base","frgion", "ion_type", "ion_type", "mod_seq", "seq", "gene"]
TYPE2LEVEL = dict(zip(TYPES, LEVELS))


def aggregate_node_properties(node, only_use_mainclust, use_fewpeps_per_protein):
    """Goes through the children and summarizes their properties to the node

    Args:
        node ([type]): [description]
        only_use_mainclust (bool, optional): [description]. Defaults to True.
    """

    if only_use_mainclust:
        childs = [x for x in node.children if x.is_included & (x.cluster ==0)]
    else:
        childs = [x for x in node.children if x.is_included]

    childs_fewpepfiltered = get_selected_nodes_for_zvalcalc(childs, use_fewpeps_per_protein, node)


    zvals = get_feature_numpy_array_from_nodes(nodes=childs_fewpepfiltered, feature_name="z_val")
    fcs =  get_feature_numpy_array_from_nodes(nodes=childs, feature_name="fc")
    cvs = get_feature_numpy_array_from_nodes(nodes=childs, feature_name="cv")
    min_intensities = get_feature_numpy_array_from_nodes(nodes = childs, feature_name = "min_intensity")
    total_intensities = get_feature_numpy_array_from_nodes(nodes = childs, feature_name = "total_intensity")
    
    min_intensity = np.median(min_intensities)
    total_intensity = np.sum(total_intensities)
    min_reps_childs = get_feature_numpy_array_from_nodes(nodes = childs, feature_name = "min_reps")
    min_reps = np.median(min_reps_childs)
    if np.isnan(min_intensity) or np.isnan(min_reps):
        Exception("values could not be determined!")

    fraction_consistent = sum([x.fraction_consistent/len(node.children) for x in childs if x.cluster ==0])



    z_sum = sum(zvals)
    p_z = NormalDist(mu = 0, sigma = np.sqrt(len(zvals))).cdf(z_sum)
    p_z = set_bounds_for_p_if_too_extreme(p_z)
    z_normed = NormalDist(mu = 0, sigma=1).inv_cdf(p_z)

    p_val = 2.0 * (1.0 - NormalDist(mu = 0, sigma = np.sqrt(len(zvals))).cdf(abs(z_sum)))
    p_val = set_bounds_for_p_if_too_extreme(p_val)

    node.z_val = z_normed
    node.p_val = p_val

    # if node.type == "frgion":
    #     node.fc = calc_weighted_fold_change_from_included_leaves_fcs(node)
    # else:
    node.fc = np.median(fcs)
    #calc_fold_change_from_included_leaves_fcs(node) ##  #np.median(fcs)#
    node.fraction_consistent = fraction_consistent
    node.cv = min(cvs)
    node.min_intensity = min_intensity
    node.total_intensity = total_intensity
    node.min_reps = min_reps

    if hasattr(node.children[0], 'predscore'):
        predscores = get_feature_numpy_array_from_nodes(nodes = childs, feature_name = "predscore")
        node.predscore = select_predscore_with_minimum_absval(predscores)
        node.cutoff = childs[0].cutoff
        node.ml_excluded = bool(abs(node.predscore)> node.cutoff)


def get_feature_numpy_array_from_nodes(nodes, feature_name ,dtype = 'float'):
    generator = (x.__dict__.get(feature_name) for x in nodes)
    return np.fromiter(generator, dtype=dtype)

def get_selected_nodes_for_zvalcalc(childs, use_fewpeps_per_protein, node):
    if use_fewpeps_per_protein and node.type == "gene":
        return filter_fewpeps_per_protein(childs)
    else:
        return childs


def filter_fewpeps_per_protein(peptide_nodes):
    peps_filtered = []
    pepnode2pval2numleaves = []
    for pepnode in peptide_nodes:
        pepleaves = [x for x in pepnode.leaves if "seq" in getattr(x,"inclusion_levels", [])]
        pepnode2pval2numleaves.append((pepnode, pepnode.z_val,len(pepleaves)))
    pepnode2pval2numleaves = sorted(pepnode2pval2numleaves, key=lambda x : abs(x[1])) #sort with lowest absolute z-val (least significant) first

    return get_median_peptides(pepnode2pval2numleaves)



import math
def get_median_peptides(pepnode2pval2numleaves):
    median_idx = math.floor(len(pepnode2pval2numleaves)/2)
    if len(pepnode2pval2numleaves)<3:
        return [x[0] for x in pepnode2pval2numleaves]
    else:
        return [x[0] for x in pepnode2pval2numleaves[:median_idx+1]]

def set_bounds_for_p_if_too_extreme(p_val):
    if p_val <aqvariables.MIN_PVAL:
        return aqvariables.MIN_PVAL
    elif p_val > 1-(aqvariables.MIN_PVAL):
        return 1- (aqvariables.MIN_PVAL)
    else:
        return p_val
    
def calc_fold_change_from_included_leaves_fcs(node):
    included_leaves = obtain_all_included_leaves(node)
    list_of_fcs = [x.fcs for x in included_leaves]
    merged_fcs = np.concatenate(list_of_fcs)
    return np.median(merged_fcs)

def calc_weighted_fold_change_from_included_leaves_fcs(node):
    included_leaves = obtain_all_included_leaves(node)
    list_of_fcs = [x.fcs for x in included_leaves]
    weights = [get_weight_of_leaf(x) for x in included_leaves]
    weighted_median = calculate_weighted_median(weights, list_of_fcs)
    return weighted_median

def get_weight_of_leaf(leaf):
    if hasattr(leaf, "predscore_fragion"):
        return 2**-leaf.predscore_fragion
    else:
        return 1

def calculate_weighted_median(weights, fcs):
    weighted_fcs = [(fc, weight) for weight, fc_list in zip(weights, fcs) for fc in fc_list]
    sorted_weighted_fcs = sorted(weighted_fcs, key=lambda x: x[0])
    sorted_fcs, sorted_weights = zip(*sorted_weighted_fcs)
    cumulative_weights = np.cumsum(sorted_weights)
    total_weight = cumulative_weights[-1]
    median_cutoff = total_weight / 2
    median_idx = np.where(cumulative_weights >= median_cutoff)[0][0]
    weighted_median = sorted_fcs[median_idx]
    return weighted_median

def obtain_all_included_leaves(node):
    list_of_included_leaves = []
    traverse_and_add_included_leaves(node, list_of_included_leaves)
    return list_of_included_leaves

def traverse_and_add_included_leaves(node, list_of_included_leaves, is_root=True):
    """
    Recursively searches for leaves from the given node, where each node in the
    path to the leaf has the 'is_included' attribute set to True, except for the initial node.
    Fills up the list_of_included_leaves with the included leaves.

    Parameters:
    node (anytree.Node): The node to start the search from.
    list_of_included_leaves (list): The list to store the included leaves in.
    is_root (bool): Indicates if the current node is the root node of the traversal.
    """

    if len(node.children) == 0:  # if the node is a leaf
        if is_root or (node.is_included and node.cluster == 0):
            list_of_included_leaves.append(node)
        return
    
    # If it's the root node or if the current node is included, then proceed to its children
    if is_root or (node.is_included and node.cluster == 0):
        for child in node.children:
            # Recursive call with is_root set to False, as we are now dealing with child nodes
            traverse_and_add_included_leaves(child, list_of_included_leaves, is_root=False)

def select_predscore_with_minimum_absval(predscores):
    abs_predscores = [abs(x) for x in predscores]
    min_value = min(abs_predscores)
    min_index = abs_predscores.index(min_value)
    return predscores[min_index]


def get_grouped_mainclust_leafs(child_nodes):
    grouped_leafs = []
    for child in child_nodes:
        child_leaves_mainclust = []
        types_previous_level = {x.type for x in child.children}
        for leafnode in child.leaves:#go through the leafs of each child
            if hasattr(leafnode, 'inclusion_levels') and not (leafnode.inclusion_levels[-1] in types_previous_level):
                continue
            child_leaves_mainclust.append(leafnode)
        if len(child_leaves_mainclust)>0:
            grouped_leafs.append(child_leaves_mainclust)
    return grouped_leafs

def select_highid_lowcv_leafs(grouped_leafs):
    grouped_leafs_lowcv = []
    for leafs in grouped_leafs:
        top_quantile_idx = math.ceil(len(leafs) * 0.2)
        leafs_repsorted = sorted(leafs, key = lambda x : x.min_reps)[:top_quantile_idx]
        leafs_repsorted_cvsorted = sorted(leafs_repsorted, key = lambda x : x.cv)
        grouped_leafs_lowcv.append([leafs_repsorted_cvsorted[0]])
    return grouped_leafs_lowcv

def select_median_fc_leafs(grouped_leafs):
    grouped_leafs_medianfc = []
    for leafs in grouped_leafs:
        leafs_fcsorted = sorted(leafs, key = lambda x : x.fc)
        grouped_leafs_medianfc.append([leafs_fcsorted[int(len(leafs_fcsorted)/2)]])
    return grouped_leafs_medianfc

def map_grouped_leafs_to_diffions(grouped_leafs, ionname2diffion):
    grouped_diffions = []
    for leafs in grouped_leafs:
        diffions = [ionname2diffion.get(x.name) for x in leafs]
        grouped_diffions.append(diffions)
    return grouped_diffions


def annotate_mainclust_leaves(childnode2clust):
    #annotate each leaf that has reached the current level with the level name, allows to visualize how the leafs are propagated
    for child in childnode2clust.keys():
        if childnode2clust.get(child)!=0:
            continue
        types_previous_level = {x.type for x in child.children}
        for leafnode in child.leaves:#annotate the leaves of each node, if they were included at this level
            if hasattr(leafnode, 'inclusion_levels'):

                if leafnode.inclusion_levels[-1] in types_previous_level: #only add a level if the previous level has also been included
                    leafnode.inclusion_levels.append(child.type)
            else:
                leafnode.inclusion_levels = [child.type]

def assign_cluster_number(type_node, childnode2clust):
    for node in type_node.children:
        if not node.is_included:
            continue
        clustid =  childnode2clust.get(node)
        node.cluster = clustid


def assign_clusterstats_to_type_node(type_node, childnode2clust):
    clust_nums = list(childnode2clust.values())
    type_node.num_clusters = len(set(clust_nums))
    type_node.num_mainclusts = sum([x==0 for x in clust_nums])
    type_node.frac_mainclust = type_node.num_mainclusts/len(clust_nums)


import scipy.stats
def assign_properties_to_base_ions(root_node, name2diffion, normed_c1, normed_c2):
    for leaf in root_node.leaves:
        log2intensities_c1 = normed_c1.ion2nonNanvals.get(leaf.name)
        log2intensities_c2 = normed_c2.ion2nonNanvals.get(leaf.name)
        diffion = name2diffion.get(leaf.name)
        leaf.fc = diffion.fc
        leaf.z_val = diffion.z_val
        leaf.fcs = get_fcs_of_leaf(log2intensities_c1, log2intensities_c2)
        leaf.fraction_consistent = 1
        original_intensities_c1 = 2**(log2intensities_c1)
        original_intensities_c2 = 2**(log2intensities_c2)
        cv_c1 = scipy.stats.variation(original_intensities_c1)
        cv_c2 = scipy.stats.variation(original_intensities_c2)
        leaf.cv = np.mean([cv_c1, cv_c2])
        leaf.min_intensity = min(sum(original_intensities_c1)/len(original_intensities_c1), sum(original_intensities_c2)/len(original_intensities_c2))
        leaf.total_intensity = get_total_intensity_if_ms2_ion(leaf, original_intensities_c1, original_intensities_c2)
        leaf.min_reps = min(len(log2intensities_c1), len(log2intensities_c2))


def get_total_intensity_if_ms2_ion(leaf, original_intensities_c1, original_intensities_c2):
    if leaf.parent.type == "ms1_isotopes":
        return 0
    else:
        return np.mean([sum(original_intensities_c1)/len(original_intensities_c1), sum(original_intensities_c2)/len(original_intensities_c2)])


def get_fcs_of_leaf(log2intensities_c1, log2intensity_c2):
    log2intensities_c1 = downsample_intensities_if_necessary(log2intensities_c1, 20)
    log2intensity_c2 = downsample_intensities_if_necessary(log2intensity_c2, 20)
    log2intensities_c1_reshaped = log2intensities_c1.reshape(-1, 1)
    fold_changes = (log2intensities_c1_reshaped - log2intensity_c2).flatten()
    return fold_changes
    

def downsample_intensities_if_necessary(intensities, max_num):
    if len(intensities) > max_num:
        rng = np.random.default_rng(42)  # Create a new random generator with the specified seed
        return rng.choice(intensities, max_num, replace=False)
    else:
        return intensities

    


def exchange_cluster_idxs(fclust_output_array):
    """The fcluster output assigns cluster numbers to the clustered elems, e.g. [1,2,1,2,2,2].
    This function here ensures that the numbers follow size of the cluster, e.g. [1,0,1,0,0,0]"""
    clustnum2count = {}
    for clustnum in fclust_output_array:
        clustnum2count[clustnum] = clustnum2count.get(clustnum, 0)+1
    clustnums = list(clustnum2count.keys())
    clustnums.sort(key = lambda x : clustnum2count.get(x), reverse= True)
    clustnum_old2clustnum_new = {clustnums[idx]: idx for idx in range(len(clustnums))}
    return [clustnum_old2clustnum_new.get(clustnum) for clustnum in fclust_output_array]


def get_fcs_ions(diffions):
    fcs = np.ones(len(diffions))
    for idx in range(len(diffions)):
        fcs[idx] = np.nanmedian([ion.fc for ion in diffions[idx]])
    return fcs



import anytree
from anytree.exporter import JsonExporter
import alphaquant.utils.utils as aqutils

from numpy import int64
from anytree import Node, iterators

def export_condpairtree_to_json(condpair_node,  results_dir):
    """exports all base roots for a given condition pair to a json file"""
    condpairname = aqutils.get_condpairname(condpair_node.name)
    results_file = f"{results_dir}/{condpairname}.iontrees.json"

    j_exporter = JsonExporter(indent=2, sort_keys=True)
    filehandle = open(results_file, "w")
    j_exporter.write(condpair_node, filehandle)
    filehandle.close()


def get_condpair_node(list_of_protein_nodes, condpair):
    condpair_node = anytree.Node(condpair) #set the condpair as node and export the whole condpair as one tree
    condpair_node.type = "condpair"
    condpair_node.level = "condpair"
    for root in list_of_protein_nodes:
        root.parent = condpair_node
    remove_unnecessary_attributes(condpair_node, ["fcs"])
    return condpair_node


def remove_unnecessary_attributes(node, attributes_to_remove):
    for attr in attributes_to_remove:
        if hasattr(node, attr):
            delattr(node, attr)

    for child in node.children:
        remove_unnecessary_attributes(child, attributes_to_remove)



import os

def get_nodes_of_type(cond1, cond2, results_folder, node_type = 'mod_seq_charge'):

    tree_sn = aqutils.read_condpair_tree(cond1, cond2, results_folder=results_folder)
    tree_sn.type = "asd"
    return anytree.findall(tree_sn, filter_= lambda x : (x.type == node_type))



def get_levelnodes_from_nodeslist(nodeslist, level):
    levelnodes = []
    for node in nodeslist:
        precursors = anytree.findall(node, filter_= lambda x : (x.type == level))
        levelnodes.extend(precursors)
    return levelnodes


def find_node_parent_at_level(node, level):
    if node.type == level:
        return node
    while node.parent is not None:
        node = node.parent
        if node.type == level:
            return node

# Cell

def check_if_node_is_included(node):
    while node.type != "gene":
        if node.cluster != 0:
            return False
        node = node.parent

    return True

def shorten_root_to_level(root, parent_level):
    for node in anytree.PreOrderIter(root):
        if node.level == parent_level:
            for child in node.children:
                child.children = tuple()
    return root



def get_parent2children_dict(tree, parent_level):
    parent2children = {}
    parent_nodes = anytree.search.findall(tree, filter_=lambda node:  node.level == parent_level)
    for parent_node in parent_nodes:
        parent2children[parent_node.name] = [child.name for child in parent_node.children]
    return parent2children

def get_parent2leaves_dict(protein):
    """Returns a dict that maps the parent node name to the names of the leaves of the parent node
    """
    parent2children = collections.defaultdict(list)
    for leave in protein.leaves:
        parent2children[leave.parent.name].append(leave.name)
    
    return dict(parent2children)

def find_max_depth( node, depth=0):
    if not node.children:
        return depth
    return max(find_max_depth(child, depth+1) for child in node.children)



def add_level_name_to_root(anynode):
    anynode.level = TYPE2LEVEL[anynode.type]
    for child in anynode.children:
        add_level_name_to_root(child)



def clone_tree(node):
    attrs = {k: v for k, v in node.__dict__.items() if not k.startswith("_")}
    
    cloned_node = anytree.Node(**attrs)
    
    for child in node.children:
        cloned_child = clone_tree(child)
        cloned_child.parent = cloned_node
    
    return cloned_node


def get_sorted_peptides_by_position_in_protein_seq(protein_node, protein_sequence):
    peptides = protein_node.children
    return sorted(peptides, key=lambda x: get_sequence_position(protein_sequence, aqutils.cut_trailing_parts_seqstring(x.name_reduced)))



def get_sequence_position(protein_seq, peptide_seq):
    return protein_seq.find(peptide_seq)


def get_sorted_peptides_by_cluster(protein_node):
    sorted_by_name = sorted(protein_node.children, key=lambda x: x.name)
    return sorted(sorted_by_name, key=lambda x: x.cluster)

def get_sorted_peptides_by_name(protein_node):
    return sorted(protein_node.children, key=lambda x: x.name)

def iterate_through_tree_levels_bottom_to_top(root_node):
    nodes_by_level = list(LevelOrderGroupIter(root_node))
    for level_nodes in reversed(nodes_by_level):
        yield level_nodes