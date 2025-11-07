import alphaquant.cluster.cluster_utils as aqcluster_utils
import anytree
import numpy as np

def apply_peptide_outlier_filtering(protnodes: list[anytree.Node]):
    regulation_score = calculate_regulation_score(protnodes)
    for protnode in protnodes:
        _determine_and_annotate_outlier_status_of_peptides(protnode, regulation_score)

        aqcluster_utils.aggregate_node_properties(protnode, only_use_mainclust=True, peptide_outlier_filtering=True)




def calculate_regulation_score(protnodes: list[anytree.Node]):
    p_values = aqcluster_utils.get_feature_numpy_array_from_nodes(nodes=protnodes, feature_name="p_val")
    abs_log2fc = np.abs(aqcluster_utils.get_feature_numpy_array_from_nodes(nodes=protnodes, feature_name="fc"))

    nonsig_mask = p_values >= 0.05
    sig_mask_005 = p_values < 0.05
    num_sig = np.sum(sig_mask_005)
    num_insig = np.sum(nonsig_mask)
    fraction_sig = num_sig / (num_sig + num_insig)

    log2fc_ratio_sig_vs_insig = np.median(abs_log2fc[sig_mask_005]) / (np.median(abs_log2fc[nonsig_mask]) + 1e-6)
    regulation_score = min(1, log2fc_ratio_sig_vs_insig * fraction_sig/10) #merges the regulation strength and the fraction of significant proteins into one score divided by to normalize it, the normalization factor corresponds to a very stongly regulated dataset
    return regulation_score




def _determine_and_annotate_outlier_status_of_peptides(protnode, regulation_score, only_use_mainclust=True):
    """
    We look at the distribution of p-values of the peptides of a protein and focus on a particular class of proteins that are potentially dominated by outliers. For these protein 1) the majority of peptides is not significant 2) there are a few peptides are more strongly significant than the majority.
    Depending on the the overall context of the experiment which is quantified by the regulation score (low regulation score means few weakly regulated proteins, high regulation score means many strongly regulated proteins), we are more or less tolerant to the outliers.

    Args:
        protnode: Protein node with peptide children
        regulation_score: Float between 0 and 1 representing overall regulation context

    Returns:
        Filtered list of inlier peptide nodes
    """
    peptide_nodes = protnode.children
    if only_use_mainclust:
        peptide_nodes = [node for node in peptide_nodes if node.cluster == 0]

    p_values = [node.p_val for node in peptide_nodes]
    median_p_val = np.median(p_values)

    if len(peptide_nodes) < 4:
        # Annotate all as inliers when too few peptides
        _annotate_peptides(peptide_nodes, is_outlier=False)
        return

    # Only check for outliers if median p-value is 0.05 or higher
    if median_p_val < 0.05:
        # Annotate all as inliers when median is already significant
        _annotate_peptides(peptide_nodes, is_outlier=False)
        return

    # Calculate threshold p-value scaled by regulation score
    threshold_p_val = median_p_val* regulation_score

    # Use list comprehensions to separate and annotate
    inlier_peptides = [node for node in peptide_nodes if node.p_val > threshold_p_val]
    outlier_peptides = [node for node in peptide_nodes if node.p_val <= threshold_p_val]

    _annotate_peptides(inlier_peptides, is_outlier=False)
    _annotate_peptides(outlier_peptides, is_outlier=True)


    return


def _annotate_peptides(peptides: list[anytree.Node], is_outlier: bool):
    """
    Annotate a list of peptides with their outlier status.

    Args:
        peptides: List of peptide nodes to annotate
        is_outlier: True if peptides are outliers, False if inliers
    """
    for peptide in peptides:
        peptide.is_outlier_peptide = is_outlier
