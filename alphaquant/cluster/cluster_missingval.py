import alphaquant.cluster.cluster_ions as aq_cluster_ions
import anytree
import alphaquant.cluster.cluster_utils as aq_cluster_utils
import numpy as np
from scipy.stats import chi2_contingency
import numpy as np
import statistics

NREP_C1 = np.nan
NREP_C2 = np.nan

def create_protnode_from_missingval_ions(gene_name, diffions, normed_c1, normed_c2, nrep_c1, nrep_c2):
    #nrep_c1 and nrep_c2 are the number of replicates in the conditions in general, not the minimum required

    global NREP_C1
    global NREP_C2
    NREP_C1 = nrep_c1
    NREP_C2 = nrep_c2

    root_node = aq_cluster_ions.create_hierarchical_ion_grouping(gene_name, diffions)
    aq_cluster_ions.add_reduced_names_to_root(root_node)


    assign_properties_to_missingval_base_ions(root_node, normed_c1, normed_c2)
    nodes_to_test = get_nodes_to_test(root_node)
    levelname_nodes_to_test = nodes_to_test[0].level
    propagate_properties_to_nodes_to_test(nodes_to_test)
    propagate_properties_from_nodes_to_test_to_root(root_node, levelname_nodes_to_test)

    return root_node

def assign_properties_to_missingval_base_ions(root_node, normed_c1, normed_c2):
    all_intensities_c1 = np.concatenate(list(normed_c1.ion2nonNanvals.values()))
    all_intensities_c2 = np.concatenate(list(normed_c2.ion2nonNanvals.values()))

    total_intensity = np.mean(np.concatenate([all_intensities_c1, all_intensities_c2]))
    lower_intensity_c1 = np.quantile(all_intensities_c1, 0.05)
    lower_intensity_c2 = np.quantile(all_intensities_c2, 0.05)

    for leaf in root_node.leaves:
        log2intensities_c1 = normed_c1.ion2nonNanvals.get(leaf.name)
        log2intensities_c2 = normed_c2.ion2nonNanvals.get(leaf.name)
        leaf.numvals_c1 = len(log2intensities_c1)
        leaf.numvals_c2 = len(log2intensities_c2)

        if leaf.numvals_c1 > leaf.numvals_c2:
            leaf.fc = np.mean(log2intensities_c1) - lower_intensity_c2
        elif leaf.numvals_c1 < leaf.numvals_c2:
            leaf.fc = lower_intensity_c1 - np.mean(log2intensities_c2)
        else:
            raise ValueError("The number of values in the conditions are equal, which should not be the case.")
        
        leaf.missingval = True
        leaf.total_intensity = total_intensity
        leaf.fraction_consistent = np.nan


def get_nodes_to_test(root_node): #get the nodes in the lowest level that is relevant for the binomial test
    if root_node.depth == 2:
        return root_node.children
    else:
        return anytree.search.findall(root_node, filter_=lambda node: node.type == "mod_seq_charge")


def propagate_properties_to_nodes_to_test(nodes_to_test): #goes through each node to test and merges the properties from it's base to the node itself
    for node in nodes_to_test:
        for level_nodes in aq_cluster_utils.iterate_through_tree_levels_bottom_to_top(node):
            if level_nodes[0].level == "base":
                continue

            for level_node in level_nodes:
                aggregate_node_properties_missingval(level_node)
                level_node.missingval = True




def propagate_properties_from_nodes_to_test_to_root(root_node, levelname_nodes_to_test):
    level_above_nodes_to_test = False
    for level_nodes in aq_cluster_utils.iterate_through_tree_levels_bottom_to_top(root_node):
        if level_nodes[0].level == levelname_nodes_to_test:
            assign_missingvals_prob_per_node(level_nodes)
            level_above_nodes_to_test = True
            continue
        if level_above_nodes_to_test:
            for level_node in level_nodes:
                aggregate_node_properties_missingval(level_node)


def aggregate_node_properties_missingval(node):
    childs = node.children
    node.numvals_c1 = np.mean([child.numvals_c1 for child in childs])
    node.numvals_c2 = np.mean([child.numvals_c2 for child in childs])
    node.fc = np.mean([child.fc for child in childs])
    node.missingval = True
    node.fraction_consistent = np.nan
    node.total_intensity = np.sum([child.total_intensity for child in childs])
    if hasattr(childs[0], "z_val"):
        node.z_val = aq_cluster_utils.sum_and_re_scale_zvalues([child.z_val for child in childs])
        node.p_val = aq_cluster_utils.transform_znormed_to_pval(node.z_val)




def assign_missingvals_prob_per_node(nodes_to_test):
    for node in nodes_to_test:
        node.p_val = perform_chisquared_test(node)
        flipped_pval = 1-0.5*node.p_val #the flipped pval is always larger than 0.5 and the closer to 1 is gets, the closer it goes to 0.5, while the smaller it gets, the closer it goes to 1. When we express this with the standard normal distribution, we are always on the right side of the distribution, so we can use the inv_cdf function to get a positive z-value equivalent to the p-value
        node.z_val = abs(statistics.NormalDist().inv_cdf(flipped_pval))
        #the p-value can be obtained again by applying the transformation: statistics.NormalDist().cdf(z)*2 - 1

def perform_chisquared_test(node_to_test):
    numvals_c1 = node_to_test.numvals_c1
    numvals_c2 = node_to_test.numvals_c2

    num_missing_c1 = NREP_C1 - numvals_c1
    num_missing_c2 = NREP_C2 - numvals_c2

    contingency_table = np.array([[numvals_c1, num_missing_c1],
                                [numvals_c2, num_missing_c2]])

    # Perform the Chi-squared test
    chi2, p, dof, expected = chi2_contingency(contingency_table)

    return p

