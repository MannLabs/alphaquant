import numpy as np
import alphaquant.cluster.cluster_utils as aq_clust_utils
import alphaquant.diffquant.doublediff_analysis as aq_diff_double

def add_proteoform_statistics_to_nodes(node2cluster : dict, take_median_ions : bool, normed_c1, normed_c2, ion2diffDist, p2z, deedpair2doublediffdist):
	"""
	Adds proteoform statistics to clustered nodes by comparing nodes within each cluster to those in the reference cluster (cluster 0). This function computes an overall p-value and fold change of fold change (fcfc) 
	for the proteoform distribution within each cluster by comparing it against the reference cluster (cluster 0). These statistics are then added to each node as attributes.
	
	The function leverages the `calc_doublediff_score` method from the `alphaquant.diffquant.doublediff_analysis` module, which requires several inputs related to ion distributions.

	Args:
		node2cluster (dict): A dictionary mapping each node to its corresponding cluster.
		take_median_ions (bool): If `True`, uses a subset of the ions (typically the median ions) for each peptide. If `False`, all ions are used.
		normed_c1: The normalized ion distribution for condition 1 (type as per `calc_doublediff_score`).
		normed_c2: The normalized ion distribution for condition 2 (type as per `calc_doublediff_score`).
		ion2diffDist: Mapping of ions to their differential distributions between the two conditions (type as per `calc_doublediff_score`).
		p2z: Mapping of p-values to z-scores (type as per `calc_doublediff_score`).
		deedpair2doublediffdist (dict[tuple[EmpiricalBackground, EmpiricalBackground]: EmpiricalBackground]): 
			A dictionary mapping pairs of empirical background distributions to the subtracted empirical background, representing double differential noise.

	Returns:
		None: The function modifies the nodes in place, annotating them with the computed proteoform statistics.
	"""

	
	if not _nodes_are_peptide_level(node2cluster):
		return
	

	cluster2nodes = _get_cluster2nodes(node2cluster)
	cluster2ions = _get_cluster2ions(cluster2nodes, take_median_ions)
	non_zero_clusters = [cluster for cluster in cluster2ions.keys() if cluster >0]
	cluster_0_ions = cluster2ions[0]
	cluster_0_nodes = cluster2nodes[0]

	for nz_cluster in non_zero_clusters:
		fcfc, pval = aq_diff_double.calc_doublediff_score(ions1=cluster_0_ions, ions2=cluster2ions[nz_cluster], 
									   normed_c1=normed_c1, normed_c2=normed_c2, ion2diffDist=ion2diffDist, p2z=p2z, 
									   deedpair2doublediffdist=deedpair2doublediffdist)
		nodes = cluster2nodes[nz_cluster]
		_annotate_nodes_with_proteoform_stats(nodes, -fcfc, pval) #fcfc direction is inverted to align with the previously defined direction of the fcfc
	
	_annotate_nodes_with_proteoform_stats(cluster_0_nodes, np.nan, np.nan)

def _nodes_are_peptide_level(node2cluster : dict):
	return all([node.type == "seq" for node in node2cluster.keys()])

def _get_cluster2nodes(node2cluster):
	cluster2nodes = {}
	for node, cluster in node2cluster.items():
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

def _annotate_nodes_with_proteoform_stats(nodes, fcfc, pval):
	for node in nodes:
		node.proteoform_fcfc = fcfc
		node.proteoform_pval = pval
