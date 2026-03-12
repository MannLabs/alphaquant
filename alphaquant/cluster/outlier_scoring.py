
import alphaquant.cluster.cluster_utils as aqclustutils
import alphaquant.utils.utils as aqutils
import numpy as np
import copy
import anytree

import alphaquant.config.config as aqconfig
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

class OutlierHandler():
    def __init__(self, condpair_tree):
        self._protnodes = condpair_tree.children

    def get_diffclust_overview_list(self):
        """_summary_

        Returns:
            ClusterDiffInfo: object containing the relevant information about two differing clusters (fcfc, peptides),
            can "reduce" protein nodes to the cluster-relevant peptides
        """
        diffclusts = []
        counter = 0
        for protnode in self._protnodes:
            counter+=1
            cluster_checker = ProtnodeClusterChecker(protnode)
            diffclusts += cluster_checker.get_diffclusts()
        return diffclusts


class ProtnodeClusterChecker():
    def __init__(self, protnode):
        self._protnode = protnode
        self._num_clusters = protnode.num_clusters

    def get_diffclusts(self):
        if not self.__check_if_multiple_clusters__():
            return []
        return self.__get_clusterdiff_info_for_each_cluster()

    def __get_clusterdiff_info_for_each_cluster(self):
        protnodes = []
        mainclust_info= self.__get_cluster_info(clustnum = 0)
        for clustnum in range(1, self._num_clusters):
            outlier_info = self.__get_cluster_info(clustnum)
            protnodes.append(self.__get_clusterdiff_info__(outlier_info, mainclust_info))
        return protnodes

    def __get_cluster_info(self, clustnum):
        mainclust_peptides = self.__get_peptides_of_cluster__(clustnum)
        return ClusterInfo(protein_name=self._protnode.name,peptide_nodes = mainclust_peptides)

    def __get_clusterdiff_info__(self, outlier_info, mainclust_info):
        return ClusterDiffInfo(mainclust_info, outlier_info)

    def __get_peptides_of_cluster__(self, clustnum):
        return [x for x in self._protnode.children if x.cluster == clustnum]

    def __check_if_multiple_clusters__(protein_node):
        return protein_node._num_clusters >1


class ClusterInfo():
    def __init__(self, protein_name,peptide_nodes):
        self.protein_name = protein_name
        self.cluster_number = list({x.cluster for x in peptide_nodes})[0]
        self.peptide_names = [x.name for x in peptide_nodes]
        self.median_fc = np.median(np.array([x.fc for x in peptide_nodes]))
        self.quality_score = self._get_quality_score(peptide_nodes)

    @staticmethod
    def _get_quality_score(peptide_nodes):
        if hasattr(peptide_nodes[0], 'ml_score'):
            return min([abs(x.ml_score) for x in peptide_nodes])
        else:
            return min(1/x.fraction_consistent for x in peptide_nodes)

class ClusterDiffInfo():
    def __init__(self, mainclust_info, outlier_info):
        self.protein_name = mainclust_info.protein_name
        self.clusterpair_id = f"{mainclust_info.cluster_number}_{outlier_info.cluster_number}"
        self.name = f"{self.protein_name}_{self.clusterpair_id}"
        self.fcdiff = abs(mainclust_info.median_fc - outlier_info.median_fc)
        self.quality_score = max(mainclust_info.quality_score, outlier_info.quality_score)
        self.outlier_peptide_names = outlier_info.peptide_names
        self.mainclust_peptide_names = mainclust_info.peptide_names
        self.peptide_names = self.mainclust_peptide_names + self.outlier_peptide_names

    def get_clusterdiff_protnode(self, protnode):
        protnode_clusterdiff = copy.deepcopy(protnode)
        self.__remove_peptides_not_in_cluster__(protnode_clusterdiff)
        self.__add_diffinfos__(protnode_clusterdiff)
        return protnode_clusterdiff

    def get_num_mainclust_peptides(self):
        return len(self.mainclust_peptide_names)

    def get_num_outlierclust_peptides(self):
        return len(self.outlier_peptide_names)

    def __remove_peptides_not_in_cluster__(self, protnode_clusterdiff):
        for peptide_node in protnode_clusterdiff.children:
            self.__remove_peptide_if_necessary__(peptide_node)

    def __add_diffinfos__(self, protnode):
        protnode.fcdiff = self.fcdiff
        protnode.quality_score = self.quality_score
        protnode.peptide_names = self.peptide_names

    def __remove_peptide_if_necessary__(self, peptide_node):
        if peptide_node.name not in self.peptide_names:
            peptide_node.parent = None



