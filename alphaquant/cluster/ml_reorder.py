import alphaquant.cluster.cluster_utils as aqcluster_utils
import numpy as np
def update_nodes_w_ml_score(protnodes):
    for prot in protnodes:
        re_order_depending_on_ml_score(prot)


def re_order_depending_on_ml_score(protnode):
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
                    clust2newclust = get_clust2newclust(child_nodes)
                    re_assign_proteoform_stats(child_nodes, clust2newclust)
                    re_order_clusters_by_ml_score(child_nodes, clust2newclust)
                    aqcluster_utils.aggregate_node_properties(type_node,only_use_mainclust=True, use_fewpeps_per_protein=True)


def get_clust2newclust(nodes):
    cluster2scores = {}
    for node in nodes:
        cluster2scores[node.cluster] = cluster2scores.get(node.cluster, [])
        cluster2scores[node.cluster].append(abs(node.ml_score))
    clusters = list(cluster2scores.keys())
    clusters.sort(key = lambda x : 1/len(cluster2scores.get(x))) 
    clusters.sort(key = lambda x : np.nanmin(cluster2scores.get(x))) #second sort preserves the order of the first sort (see test_ml_reorder.py)
    clust2newclust = { clusters[x] :x for x in range(len(clusters))}
    return clust2newclust

def re_assign_proteoform_stats(nodes, clust2newclust):
	if nodes[0].level !="sequence":
		return
	zero_cluster_has_changed = clust2newclust[0] != 0
	if zero_cluster_has_changed:
		change_pformstats_from_old_to_new_cluster(nodes, 0, clust2newclust[0])


def change_pformstats_from_old_to_new_cluster(nodes, zero_cluster, new_zero_cluster):
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


def re_order_clusters_by_ml_score(nodes, clust2newclust):
    for node in nodes:
        node.cluster =clust2newclust.get(node.cluster)



