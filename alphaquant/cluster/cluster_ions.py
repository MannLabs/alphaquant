import scipy.spatial.distance
import scipy.cluster.hierarchy
import alphaquant.cluster.cluster_utils as aqcluster_utils
import alphaquant.cluster.cluster_sorting as aq_cluster_sorting
import alphaquant.diffquant.diffutils as aqutils
import statsmodels.stats.multitest as multitest
import numpy as np

import alphaquant.config.config as aqconfig
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)


REGEX_FRGIONS_ISOTOPES = [[("(SEQ.*MOD.*CHARGE.*FRG)(ION.*)", "frgion"), ("(SEQ.*MOD.*CHARGE.*MS1)(ISO.*)", "ms1_isotopes"), ("(SEQ.*MOD.*CHARGE.*PREC)(URSOR.*)", "precursor")], [("(SEQ.*MOD.*CHARGE.*)(_FRG.*|_MS1.*|_PREC.*)", "mod_seq_charge")], [("(SEQ.*MOD.*)(CHARGE.*)", "mod_seq")], [("(SEQ.*)(MOD.*)", "seq")]]
LEVEL_NAMES = ['ion_type', 'mod_seq_charge', 'mod_seq', 'seq']
MAPPING_DICT = {'SEQ':'seq', 'MOD':'mod_seq', 'CHARGE':'mod_seq_charge', 'MS1ISOTOPES':'ms1_isotopes','FRGION':'frgion', 'PRECURSOR' : 'precursor'}
FCDIFF_CUTOFF_CLUSTERMERGE = 0
LEVEL2PVALTHRESH = {'ion_type':0.01, 'mod_seq_charge':0.01, 'mod_seq':1e-20, 'seq':1e-20} #the pval threshold is only set at the gene level, the rest of the levels are set as specified here. The threshold applies to the children of the node





def get_scored_clusterselected_ions(gene_name, diffions, normed_c1, normed_c2, ion2diffDist, p2z, deedpair2doublediffdist, pval_threshold_basis, fcfc_threshold, take_median_ion,
                                    fcdiff_cutoff_clustermerge):
    #typefilter = TypeFilter('successive')

    global FCDIFF_CUTOFF_CLUSTERMERGE
    FCDIFF_CUTOFF_CLUSTERMERGE = fcdiff_cutoff_clustermerge

    diffions = sorted(diffions, key = lambda x : x.name)
    name2diffion = {x.name : x for x in diffions}
    root_node = create_hierarchical_ion_grouping(gene_name, diffions)
    add_reduced_names_to_root(root_node)
    #LOGGER.info(anytree.RenderTree(root_node))
    root_node_clust = cluster_along_specified_levels(root_node, name2diffion, normed_c1, normed_c2, ion2diffDist, p2z, deedpair2doublediffdist, pval_threshold_basis, fcfc_threshold, take_median_ion)

    level_sorted_nodes = [[node for node in children] for children in anytree.ZigZagGroupIter(root_node_clust)]
    level_sorted_nodes.reverse() #the base nodes are first

    root_node_lvl = level_sorted_nodes[-1]
    if len(root_node_lvl)!=1:
        Exception("there should be only one root node!")

    root_node_annot =root_node_lvl[0]
    return root_node_annot


import anytree
import re
def create_hierarchical_ion_grouping(gene_name, diffions):
    #regex patterns sorted from bottom to top in the following way list(list(tuple(pattern, name))): first instance of list represents the level of the tree, second instance represents the different nodes available on this level (for example FRgIon, MS1 are on the same level)

    nodes = [anytree.Node(x.name, type = "base", level = "base",cluster = -1, is_included = True) for x in diffions] #construct the base nodes

    for level_idx, level in enumerate(REGEX_FRGIONS_ISOTOPES):
        name2node = {}
        for pattern2name in level:
            for node in nodes:
                if (re.match(pattern2name[0], node.name)):
                    level_name = LEVEL_NAMES[level_idx]
                    m = re.match(pattern2name[0], node.name)
                    matching_name = m.group(1)
                    name2node[matching_name] = name2node.get(matching_name, anytree.Node(matching_name,  type = pattern2name[1], level = level_name,cluster = -1, is_included = True))
                    node.parent = name2node.get(matching_name)

        if len(name2node.keys())>0:
            nodes = list(name2node.values())

    root_node = anytree.Node(gene_name, type = "gene", level = "gene",cluster = 0, is_included = True)

    for node in nodes:
        node.parent = root_node

    return root_node

def add_reduced_names_to_root(node):
    for child in node.children:
        add_reduced_names_to_root(child)
    if node.parent:
        node.name_reduced = node.name.replace(node.parent.name, "")
    else:
        node.name_reduced = node.name
    

import pandas as pd
def cluster_along_specified_levels(root_node, ionname2diffion, normed_c1, normed_c2, ion2diffDist, p2z, deedpair2doublediffdist, pval_threshold_basis, fcfc_threshold, take_median_ion):#~60% of overall runtime
    #typefilter object specifies filtering and clustering of the nodes
    aqcluster_utils.assign_properties_to_base_ions(root_node, ionname2diffion, normed_c1, normed_c2)
    
    for level_nodes in  aqcluster_utils.iterate_through_tree_levels_bottom_to_top(root_node):
        nodetypes_at_level = list(set([node.type for node in level_nodes]))
        if nodetypes_at_level == ["base"]:
            continue
        for node_type in nodetypes_at_level:
            type_nodes = [x for x in level_nodes if x.type == node_type] #this gets e.g. all the precursors
            
            if len(type_nodes)==0:
                continue
            for type_node in type_nodes: #this goes through each precursor individually and clusters the children
                child_nodes = type_node.children
                grouped_mainclust_leafs = aqcluster_utils.get_grouped_mainclust_leafs(child_nodes) #leafs are excluded if they are not in the main cluster
                
                if len(grouped_mainclust_leafs)==0: #this means the leafs were previously excluded
                    exclude_node(type_node)
                    continue

                if len(grouped_mainclust_leafs)==1:
                    childnode2clust = get_childnode2clust_for_single_ion(type_node)
                else:
                    if take_median_ion:
                        grouped_mainclust_leafs = aqcluster_utils.select_median_fc_leafs(grouped_mainclust_leafs)
                    diffions = aqcluster_utils.map_grouped_leafs_to_diffions(grouped_mainclust_leafs, ionname2diffion)
                    childnode2clust = find_fold_change_clusters(type_node, diffions, normed_c1, normed_c2, ion2diffDist, p2z, deedpair2doublediffdist, pval_threshold_basis, fcfc_threshold) #the clustering is performed on the child nodes
                    childnode2clust = merge_similar_clusters_if_applicable(childnode2clust, type_node, fcdiff_cutoff_clustermerge = FCDIFF_CUTOFF_CLUSTERMERGE)
                    childnode2clust = aq_cluster_sorting.decide_cluster_order(childnode2clust)
                
                aqcluster_utils.assign_clusterstats_to_type_node(type_node, childnode2clust)
                aqcluster_utils.annotate_mainclust_leaves(childnode2clust)
                aqcluster_utils.assign_cluster_number(type_node, childnode2clust)
                aqcluster_utils.aggregate_node_properties(type_node,only_use_mainclust=True, use_fewpeps_per_protein=True)

    return root_node

def get_childnode2clust_for_single_ion(type_node):
    type_node.num_clusters = 1
    type_node.num_mainclusts = 1
    type_node.frac_mainclust = 1
    return {type_node.children[0]: 0}


def find_fold_change_clusters(type_node, diffions, normed_c1, normed_c2, ion2diffDist, p2z, deedpair2doublediffdist, pval_threshold_basis, fcfc_threshold):
    """Compares the fold changes of the ions corresponding to the nodes that are compared and returns the set of ions with consistent fold changes.

    Args:
        diffions (list[list[ionnames]]): contains the sets of ions to be tested, for example [[fragion1_precursor1, fragion2_precursor1, fragion3_precursor1],[fragion1_precursor2],[fragion1_precursor3, fragion2_precursor3]]. The ions are assumed to be similar in type (e.g. fragment, precursor)!
        normed_c1 (ConditionBackground): [description]
        normed_c2 (ConditionBackground): [description]
        ion2diffDist (dict(ion : SubtractedBackground)): [description]
        p2z ([type]): [description]
        deedpair2doublediffdist ([type]): [description]
        fc_threshold (float, optional): [description]. Defaults to 0.
        pval_threshold_basis (float, optional): the threshold at which to merge peptides at the gene level. Defaults to 0.01
    """

    pval_threshold_basis = get_pval_threshold_basis(type_node, pval_threshold_basis)
    diffions_idxs = [[x] for x in range(len(diffions))]
    diffions_fcs = aqcluster_utils.get_fcs_ions(diffions)
    #mt_corrected_pval_thresh = pval_threshold_basis/len(diffions)
    condensed_similarity_matrix = scipy.spatial.distance.pdist(diffions_idxs, lambda idx1, idx2: evaluate_similarity(idx1[0], idx2[0], diffions, diffions_fcs, normed_c1, normed_c2, ion2diffDist,p2z, 
                                                                                                   deedpair2doublediffdist, fcfc_threshold)) #gives p-values of the pairwise comparisons of the ions
    condensed_similarity_matrix_mt_corrected = get_multiple_testing_corrected_condensed_similarity_matrix(condensed_similarity_matrix)
    condensed_distance_matrix_mt_corrected = 1/condensed_similarity_matrix_mt_corrected
    
    after_clust = scipy.cluster.hierarchy.ward(condensed_distance_matrix_mt_corrected)
    clustered = scipy.cluster.hierarchy.fcluster(after_clust, 1/(pval_threshold_basis), criterion='distance')
    clustered = aqcluster_utils.exchange_cluster_idxs(clustered)

    childnode2clust = [(type_node.children[ion_idx],clust_idx) for ion_idx, clust_idx in zip(list(range(len(clustered))),clustered)]
    childnode2clust = sorted(childnode2clust, key = lambda x : x[0].name) #sort list for reproducibility

    return childnode2clust

def get_pval_threshold_basis(type_node, pval_threshold_basis): #the pval threshold is only set at the gene level, the rest of the levels are set as specified in the LEVEL2PVALTHRESH dictionary
    if type_node.level == "gene":
        return pval_threshold_basis
    else:
        return LEVEL2PVALTHRESH.get(type_node.level, 0.2)
    
def get_multiple_testing_corrected_condensed_similarity_matrix(condensed_distance_matrix: np.array):
    """
    condensed_distance_matrix contains all p-values of the pairwise comparisons of the ions. They are by definition dependent.
    
    Args:
    condensed_distance_matrix (np.array): Condensed distance matrix containing p-values of pairwise comparisons.
    
    Returns:
    np.array: Corrected condensed distance matrix.
    """
    # Apply Benjamini-Yekutieli correction
    _, corrected_pvalues, _, _ = multitest.multipletests(condensed_distance_matrix, method='fdr_by')
    
    # Return the corrected condensed matrix
    return corrected_pvalues


def merge_similar_clusters_if_applicable(childnode2clust, type_node, fcdiff_cutoff_clustermerge = 0.5):
    if type_node.level == "gene":
        return merge_similar_clusters(childnode2clust, fcdiff_cutoff_clustermerge)
    else:
        return childnode2clust


def merge_similar_clusters(childnode2clust, fcdiff_cutoff_clustermerge = 0.5):
    clust2childnodes = aqutils.invert_tuple_list_w_nonunique_values(childnode2clust)

    if len(clust2childnodes.keys())==1:
        return childnode2clust

    clust2fc = {}
    for clust, childnodes in clust2childnodes.items():
        clust2fc[clust] = np.median([x.fc for x in childnodes])

    clusters = list(clust2fc.keys())
    clust_idxs = [[x] for x in range(len(clusters))]

    condensed_distance_matrix = scipy.spatial.distance.pdist(clust_idxs, lambda idx1, idx2: compare_fcdistance(clusters, idx1, idx2, clust2fc))
    after_clust = scipy.cluster.hierarchy.complete(condensed_distance_matrix)
    clustered = scipy.cluster.hierarchy.fcluster(after_clust, fcdiff_cutoff_clustermerge, criterion='distance')

    childnode2clust = update_childnode2clust(childnode2clust, clusters, clustered)

    return childnode2clust


def compare_fcdistance(clusters, idx1, idx2, clust2fc):
    clust1 = clusters[idx1[0]]
    clust2 = clusters[idx2[0]]

    fc1 = clust2fc.get(clust1)
    fc2 = clust2fc.get(clust2)

    return abs(fc1-fc2)



def update_childnode2clust(childnode2clust, old_clusters, new_clusters):
    old2new = dict(zip(old_clusters, new_clusters))
    childnode2clust_new = []
    for childnode, old_clust in childnode2clust:
        new_clust = old2new[old_clust]
        childnode2clust_new.append((childnode, new_clust))
    return childnode2clust_new
    




# Cell
import statistics
import alphaquant.diffquant.doublediff_analysis as aqdd
import numpy as np
def evaluate_similarity(idx1, idx2, diffions, fcs, normed_c1, normed_c2, ion2diffDist, p2z, deedpair2doublediffdist, fcfc_threshold):
    """
    Evaluate the statistical similarity between two sets of ions based on their properties and fold changes.

    This function calculates a p-value representing the statistical similarity between two sets of ions,
    testing the null hypothesis that the two sets are not significantly different.

    Args:
        idx1 (int): Index of the first set of ions in the diffions list.
        idx2 (int): Index of the second set of ions in the diffions list.
        diffions (list): List of ion objects, each containing a 'name' attribute.
        fcs (list): List of fold change values corresponding to each set of ions.
        normed_c1 (array-like): Background distributions for condition 1.
        normed_c2 (array-like): Background distributions for condition 2.
        ion2diffDist (dict): Mapping of ion pairs to their difference distributions.
        p2z (function): Function to convert p-values to z-scores.
        deedpair2doublediffdist (dict): Mapping of ion pairs to their double difference distributions.
        fcfc_threshold (float): Threshold for considering fold changes as similar.

    Returns:
        float: A p-value where higher values suggest greater similarity between ion sets.
               Returns 0.99 for fold changes below fcfc_threshold.
    """
    ions1 = [x.name for x in diffions[idx1]]
    ions2 = [x.name for x in diffions[idx2]]
    fc1 = fcs[idx1]
    fc2 = fcs[idx2]

    if abs((fc1-fc2)) < fcfc_threshold:
        return 0.99 #

    fcfc, pval = aqdd.calc_doublediff_score(ions1, ions2, normed_c1, normed_c2,ion2diffDist,p2z, deedpair2doublediffdist)
    return (pval + 1e-17)



# Cell
import anytree
def exclude_node(node):
    node.is_included = False
    for descendant in node.descendants:
        descendant.is_included = False


# Cell
import numpy as np
def update_nodes_w_ml_score(protnodes):
    for prot in protnodes:
        re_order_depending_on_predscore(prot)


def re_order_depending_on_predscore(protnode):
    for level_nodes in aqcluster_utils.iterate_through_tree_levels_bottom_to_top(protnode):
        node_types = list(set([node.type for node in level_nodes])) # a certain tree level can contain different types of nodes, for example level ion_type has ms1 isotopes and frgions
        if node_types == ["base"]:
            continue
        for node_type in node_types:
            type_nodes = [x for x in level_nodes if x.type == node_type]
            if len(type_nodes)==0:
                continue
            for type_node in type_nodes: #go through the nodes, re-order the children. Propagate the values from the newly ordered children to the type node
                child_nodes = type_node.children
                had_predscore = hasattr(child_nodes[0], 'predscore')
                if had_predscore:
                    re_order_clusters_by_predscore(child_nodes)
                    aqcluster_utils.aggregate_node_properties(type_node,only_use_mainclust=True, use_fewpeps_per_protein=True)



def re_order_clusters_by_predscore(nodes):
    cluster2scores = {}
    for node in nodes:
        cluster2scores[node.cluster] = cluster2scores.get(node.cluster, [])
        cluster2scores[node.cluster].append(abs(node.predscore))
    clusters = list(cluster2scores.keys())
    clusters.sort(key = lambda x : 1/len(cluster2scores.get(x))) 
    clusters.sort(key = lambda x : np.nanmin(cluster2scores.get(x))) 
    clust2newclust = { clusters[x] :x for x in range(len(clusters))}
    for node in nodes:
        node.cluster =clust2newclust.get(node.cluster)
