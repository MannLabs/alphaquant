import alphaquant.cluster.cluster_utils as aqcluster_utils
import numpy as np
import anytree


def update_nodes_w_ml_score(protnodes : list[anytree.Node]):
    """
    Update and re-order clusters within protein nodes based on ML scores.

    This function serves as the entry point for updating a list of protein nodes.
    It iterates through the provided nodes and calls the re-ordering function
    for each one.

    Args:
        protnodes (list[anytree.Node]): A list of protein nodes to be processed.

    Returns:
        None
    """

    for prot in protnodes:
        _re_order_depending_on_ml_score(prot)


def _re_order_depending_on_ml_score(protnode : anytree.Node):
    """
    Reorder clusters in a protein node tree based on machine learning scores.

    This function traverses a tree of protein nodes from bottom to top, reorganizing
    the clusters within each node based on their machine learning (ml) scores. The
    process involves:

    1. Evaluating the ml_score of each cluster within a node.
    2. Selecting the cluster with the best ml_score as the first cluster.
    3. Updating the parent node's properties if the first cluster has changed.

    This reordering ensures that the clusters with the best ml_score are prioritized in the tree structure.

    Args:
        protnode (anytree.Node): The protein node to be processed.

    Returns:
        None
    """
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
                had_ml_score = hasattr(child_nodes[0], 'ml_score')
                if had_ml_score:
                    clust2newclust = _get_clust2newclust(child_nodes)
                    _re_assign_proteoform_stats(child_nodes, clust2newclust)
                    _re_order_clusters_by_ml_score(child_nodes, clust2newclust)
                    aqcluster_utils.aggregate_node_properties(type_node,only_use_mainclust=True, use_fewpeps_per_protein=True)


def _get_clust2newclust(nodes: list[anytree.Node]) -> dict[int, int]:
    """
    Creates a mapping from current cluster IDs to new cluster IDs based on ML scores.

    This function takes a list of nodes at the same level in a tree, each with a 
    cluster ID and an ML score. It reorders the clusters based on their ML scores 
    and size, assigning new cluster IDs.

    Args:
        nodes (list[anytree.Node]): A list of tree nodes, each representing a data point. 
                                    Each node should have 'cluster' and 'ml_score' attributes.

    Returns:
        dict[int, int]: A dictionary mapping current cluster IDs to new cluster IDs.
                        Keys are the original cluster IDs, values are the new cluster IDs.

    Note:
        - Clusters are first sorted by size (ascending), then by minimum absolute ML score.
        - The sorting preserves the order of the first sort (see test_ml_reorder.py).
        - New cluster IDs are assigned sequentially starting from 0.
    """
    cluster2scores = {}
    for node in nodes:
        cluster2scores[node.cluster] = cluster2scores.get(node.cluster, [])
        cluster2scores[node.cluster].append(abs(node.ml_score))
    clusters = list(cluster2scores.keys())
    clusters.sort(key = lambda x : 1/len(cluster2scores.get(x))) 
    clusters.sort(key = lambda x : np.nanmin(cluster2scores.get(x))) #second sort preserves the order of the first sort (see test_ml_reorder.py)
    clust2newclust = { clusters[x] :x for x in range(len(clusters))}
    return clust2newclust

def _re_assign_proteoform_stats(nodes: list[anytree.Node], clust2newclust: dict[int, int]):
    """
    Reassign proteoform statistics when the main cluster changes.

    The proteoform statistics (FDR and FCFC) of a whole cluster are always assigned 
    relative to the main cluster. This function reassigns these stats when the main 
    cluster changes.

    Args:
        nodes (list[anytree.Node]): List of nodes representing the clusters.
        clust2newclust (dict[int, int]): Mapping of old cluster IDs to new cluster IDs.

    Returns:
        None
    """
    if nodes[0].level !="sequence":
        return
    zero_cluster_has_changed = clust2newclust[0] != 0
    if zero_cluster_has_changed:
        _change_pformstats_from_old_to_new_cluster(nodes, 0, clust2newclust[0])


def _change_pformstats_from_old_to_new_cluster(nodes : list[anytree.Node], zero_cluster : int, new_zero_cluster : int):
    """
    Change proteoform statistics from the old main cluster to the new main cluster.

    This function swaps the proteoform statistics (FCFC and p-value) between the old 
    main cluster (zero cluster) and the new main cluster.

    Args:
        nodes (list[anytree.Node]): List of nodes representing the clusters.
        zero_cluster (int): ID of the old main cluster (usually 0).
        new_zero_cluster (int): ID of the new main cluster.

    Returns:
        None
    """
    nodes_zero = [node for node in nodes if node.cluster == zero_cluster]
    nodes_new_zero = [node for node in nodes if node.cluster == new_zero_cluster]
    
    proteoform_fcfc_old = nodes_zero[0].proteoform_fcfc
    proteoform_pval_old = nodes_new_zero[0].proteoform_pval
    
    proteoform_fcfc_new = nodes_new_zero[0].proteoform_fcfc
    proteoform_pval_new = nodes_new_zero[0].proteoform_pval
    
    for node in nodes_zero:
        node.proteoform_fcfc = proteoform_fcfc_new
        node.proteoform_pval = proteoform_pval_new

    for node in nodes_new_zero:
        node.proteoform_fcfc = proteoform_fcfc_old
        node.proteoform_pval = proteoform_pval_old


def _re_order_clusters_by_ml_score(nodes : list[anytree.Node], clust2newclust : dict[int, int]):
    """
    Reorder clusters based on a new ordering defined by machine learning scores.

    This function updates the cluster ID of each node according to the new ordering 
    provided in the clust2newclust dictionary.

    Args:
        nodes (list[anytree.Node]): List of nodes representing the clusters.
        clust2newclust (dict[int, int]): Mapping of old cluster IDs to new cluster IDs.

    Returns:
        None
    """
    for node in nodes:
        node.cluster =clust2newclust.get(node.cluster)



