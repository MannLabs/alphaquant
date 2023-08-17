
# Cell
import anytree
import alphaquant.diff_analysis as aqdiff
import alphaquant.diffquant_utils as aqutils
from statistics import NormalDist
import statistics
import numpy as np

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

    if use_fewpeps_per_protein and node.type == "gene":
        childs = filter_fewpeps_per_protein(childs)


    zvals = get_feature_numpy_array_from_nodes(nodes=childs, feature_name="z_val")
    fcs =  get_feature_numpy_array_from_nodes(nodes=childs, feature_name="fc")
    cvs = get_feature_numpy_array_from_nodes(nodes=childs, feature_name="cv")
    min_intensities = get_feature_numpy_array_from_nodes(nodes = childs, feature_name = "min_intensity")
    min_intensity = np.median(min_intensities)
    min_reps_childs = get_feature_numpy_array_from_nodes(nodes = childs, feature_name = "min_reps")
    min_reps = np.median(min_reps_childs)
    if np.isnan(min_intensity) or np.isnan(min_reps):
        Exception("values could not be determined!")

    fraction_consistent = sum([x.fraction_consistent/len(node.children) for x in childs if x.cluster ==0])



    z_sum = sum(zvals)
    p_z = NormalDist(mu = 0, sigma = np.sqrt(len(zvals))).cdf(z_sum)
    p_z = set_bounds_for_p_if_too_extreme(p_z)
    z_normed = NormalDist(mu = 0, sigma=1).inv_cdf(p_z)
    if z_normed <-8.3:
        Exception("not in alignment with bounded pval")
    if z_normed > 8.3:
        Exception("not in alignment with bounded pval")

    p_val = max(1e-16, 2.0 * (1.0 - NormalDist(mu = 0, sigma = np.sqrt(len(zvals))).cdf(abs(z_sum))))

    node.z_val = z_normed
    node.p_val = p_val
    node.fc = np.median(fcs)
    node.fraction_consistent = fraction_consistent
    node.cv = min(cvs)
    node.min_intensity = min_intensity
    node.min_reps = min_reps

    if hasattr(node.children[0], 'predscore'):
        predscores = [x.predscore for x in childs]
        node.predscore = select_predscore_with_minimum_absval(predscores)
        node.cutoff = childs[0].cutoff
        node.ml_excluded = bool(abs(node.predscore)> node.cutoff)

def get_feature_numpy_array_from_nodes(nodes, feature_name ,dtype = 'float'):
    generator = (x.__dict__.get(feature_name) for x in nodes)
    return np.fromiter(generator, dtype=dtype)

def filter_fewpeps_per_protein(peptide_nodes):
    peps_filtered = []
    pepnode2pval2numleaves = []
    for pepnode in peptide_nodes:
        pepleaves = [x for x in pepnode.leaves if "seq" in getattr(x,"inclusion_levels", [])]
        pepnode2pval2numleaves.append((pepnode, pepnode.p_val,len(pepleaves)))
    pepnode2pval2numleaves = sorted(pepnode2pval2numleaves, key=lambda x : x[1], reverse=True) #sort with highest p-val (least significant) first

    return get_median_peptides(pepnode2pval2numleaves)


def set_bounds_for_p_if_too_extreme(p_val):
    if p_val <1e-16:
        return 1e-16
    elif p_val > 1-(1e-16):
        return 1- (1e-16)
    else:
        return p_val

import math
def get_median_peptides(pepnode2pval2numleaves):
    median_idx = math.floor(len(pepnode2pval2numleaves)/2)
    if len(pepnode2pval2numleaves)<3:
        return [x[0] for x in pepnode2pval2numleaves]
    else:
        return [x[0] for x in pepnode2pval2numleaves[:median_idx+1]]

def select_predscore_with_minimum_absval(predscores):
    abs_predscores = [abs(x) for x in predscores]
    min_value = min(abs_predscores)
    min_index = abs_predscores.index(min_value)
    return predscores[min_index]

def get_mainclust_leaves(child_nodes, ionname2diffion):
    grouped_leafs = []
    for child in child_nodes:
        child_leaves_mainclust = []
        types_previous_level = {x.type for x in child.children}
        for leafnode in child.leaves:#go through the leafs of each child
            if hasattr(leafnode, 'inclusion_levels') and not (leafnode.inclusion_levels[-1] in types_previous_level):
                continue
            child_leaves_mainclust.append(leafnode)
        child_leafs_diffions = [ionname2diffion.get(x.name) for x in child_leaves_mainclust] #map the leaf names to the diffion objetcs
        if len(child_leafs_diffions)>0:
            grouped_leafs.append(child_leafs_diffions)
    return grouped_leafs


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
def assign_fcs_to_base_ions(root_node, name2diffion, normed_c1, normed_c2):
    for leaf in root_node.leaves:
        leaf.fc = name2diffion.get(leaf.name).fc
        leaf.z_val = name2diffion.get(leaf.name).z_val
        leaf.fraction_consistent = 1
        original_intensities_c1 = 2**(normed_c1.ion2nonNanvals.get(leaf.name))
        original_intensities_c2 = 2**(normed_c2.ion2nonNanvals.get(leaf.name))
        cv_c1 = scipy.stats.variation(original_intensities_c1)
        cv_c2 = scipy.stats.variation(original_intensities_c2)
        leaf.cv = min(cv_c1, cv_c2)
        leaf.min_intensity = min(sum(original_intensities_c1)/len(original_intensities_c1), sum(original_intensities_c2)/len(original_intensities_c2))
        leaf.min_reps = min(len(normed_c1.ion2nonNanvals.get(leaf.name)), len(normed_c2.ion2nonNanvals.get(leaf.name)) )


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


def get_diffresults_from_clust_root_node(root_node):
    pval = root_node.p_val
    fc = root_node.fc
    ions_included = [x.name for x in root_node.leaves if x.is_included]
    consistency_score = root_node.fraction_consistent * len(root_node.leaves)
    return pval, fc, consistency_score, ions_included

import anytree
from anytree.exporter import JsonExporter
import alphaquant.diffquant_utils as aqutils

from numpy import int64
from anytree import Node, iterators

def export_roots_to_json(rootlist, condpair, results_dir):
    """exports all base roots for a given condition pair to a json file"""
    condpairname = aqutils.get_condpairname(condpair)
    condpair_node = anytree.Node(condpair) #set the condpair as node and export the whole condpair as one tree
    for root in rootlist:
        root.parent = condpair_node
    results_file = f"{results_dir}/{condpairname}.iontrees.json"


    # Assuming 'root' is the root of your anytree object
    nodes = list(iterators.PreOrderIter(root))
    for node in nodes:
        for var in vars(node):
            value = getattr(node, var)
            if isinstance(value, int64):
                print(f"Node {node.name} has an {var} with value {value} of type int64")

    j_exporter = JsonExporter(indent=2, sort_keys=True)
    filehandle = open(results_file, "w")
    j_exporter.write(condpair_node, filehandle)
    filehandle.close()