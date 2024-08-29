import numpy as np
import alphaquant.cluster.cluster_utils as aq_clust_utils
import alphaquant.diffquant.doublediff_analysis as aq_diff_double

def add_proteoform_statistics_to_nodes(childnode2clust_ordered : dict, take_median_ions : bool, normed_c1, normed_c2, ion2diffDist, p2z, deedpair2doublediffdist):

	if next(iter(childnode2clust_ordered)).type != "sequence":
		return
	 
	cluster2nodes = _get_cluster2nodes(childnode2clust_ordered)
	cluster2ions = _get_cluster2ions(cluster2nodes, take_median_ions)
	non_zero_clusters = [cluster for cluster in cluster2ions.keys() if cluster >0]
	cluster_0_ions = cluster2ions[0]
	cluster_0_nodes = cluster2nodes[0]

	for nz_cluster in non_zero_clusters:
		fcfc, pval = aq_diff_double.calc_doublediff_score(ions1=cluster_0_ions, ions2=cluster2ions[nz_cluster], 
									   normed_c1=normed_c1, normed_c2=normed_c2, ion2diffDist=ion2diffDist, p2z=p2z, 
									   deedpair2doublediffdist=deedpair2doublediffdist)
		nodes = cluster2nodes[nz_cluster]
		_annotate_nodes_with_proteoform_stats(fcfc, pval, nodes)
	
	_annotate_nodes_with_proteoform_stats(np.nan, np.nan, cluster_0_nodes)
		

def _get_cluster2nodes(childnode2clust_ordered):
	cluster2nodes = {}
	for node, cluster in childnode2clust_ordered.items():
		if cluster not in cluster2nodes:
			cluster2nodes[cluster] = []
		cluster2nodes[cluster].append(node)
	return cluster2nodes

def _get_cluster2ions(cluster2nodes, take_median_ions):
	cluster2ions = {}
	for cluster, nodes in cluster2nodes.items():
		cluster2ions[cluster] = []
		for node in nodes:
			leavenames = _get_leavenames_from_node(node, take_median_ions)
			cluster2ions[cluster].extend(leavenames)
	return cluster2ions

def _get_leavenames_from_node(node, take_median_ion):
	leaves = node.leaves
	if take_median_ion:
		middle_leaves = aq_clust_utils.select_middle_leafs(leaves)
	else:
		middle_leaves = leaves

	leavenames = [leaf.name for leaf in middle_leaves]
	return leavenames

def _annotate_nodes_with_proteoform_stats(fcfc, pval, nodes):
	for node in nodes:
		node.proteoform_fcfc = fcfc
		node.proteoform_pval = pval
