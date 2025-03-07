
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



class ProtNodeRetriever():
    @staticmethod
    def get_protnodes_from_diffclust_list(condpair_tree, diffclust_list):
        proteins_wanted = {x.protein_name for x in diffclust_list}
        protnodes = list(filter(lambda x : x.name in proteins_wanted, condpair_tree.children))
        return protnodes

# Cell

class OutlierPeptideLoader():
    def __init__(self, condpair_tree):
        self._condpair_tree = condpair_tree
        self.outlier_peptides = []
        self._add_outlier_peptides()

    def _add_outlier_peptides(self):
        for protnode in self._condpair_tree.children:
            nodechecker = ProtnodeClusterCheckerPeptideInfos(protnode)
            self.outlier_peptides += nodechecker.get_outlier_peptide_infos()


class ProtnodeClusterCheckerPeptideInfos(ProtnodeClusterChecker):
    def __init__(self, protnode):
        super().__init__(protnode)
        self._outlier_peptide_infos = []

    def get_outlier_peptide_infos(self):
        diffclusts = self.get_diffclusts()
        for clusterdiffinfo in diffclusts:
            self._update_outlier_peptide_infos(clusterdiffinfo)
        return self._outlier_peptide_infos

    def _update_outlier_peptide_infos(self, clusterdiffinfo):
        peptide_nodes = self._get_outlier_peptide_nodes(clusterdiffinfo)
        for peptide_node in peptide_nodes:
            self._outlier_peptide_infos.append(OutlierPeptideInfo(peptide_node))

    def _get_outlier_peptide_nodes(self, clusterdiffinfo):
        peptide_names = set(clusterdiffinfo.outlier_peptide_names)
        return anytree.findall(self._protnode, filter_= lambda x : x.name in peptide_names, maxlevel=2)


class ProteinInfo():
    def __init__(self, peptide_node):
        self.protein_fc = self._get_protein_fc(peptide_node)

    def _get_protein_fc(self, peptide_node):
        return aqclustutils.find_node_parent_at_level(peptide_node, "gene").fc


class OutlierPeptideInfo(ProteinInfo):
    def __init__(self, peptide_node):
        super().__init__(peptide_node)
        self._peptide_node = peptide_node
        self.peptide_sequence = peptide_node.name
        self.fc = peptide_node.fc
        self.quality_score = self._get_quality_score(peptide_node)
        self.protnormed_fc = None
        self.num_mainclust_peptides = self._get_number_mainclust_peptides()
        self._calc_protnormed_fc()

    def _get_quality_score(self, peptide_node):
        has_ml_score = hasattr(peptide_node, 'ml_score')
        if has_ml_score:
            return abs(peptide_node.ml_score)
        else:
            return 1/peptide_node.fraction_consistent

    def _calc_protnormed_fc(self):
        self.protnormed_fc = self.fc - self.protein_fc

    def _get_number_mainclust_peptides(self):
        samelevel_nodes = self._peptide_node.parent.children
        mainclust_nodes = filter(lambda x : x.cluster ==0, samelevel_nodes)
        return len(list(mainclust_nodes))


# Cell
import anytree


class ModifiedPeptideLoader():
    def __init__(self, condpair_tree, specific_modification = "[Phospho (STY)]"):
        self.specific_modification = specific_modification
        self.condpair_tree = condpair_tree
        self._pepname2modpep = {}
        self._load_modified_peptides_from_tree()

    def get_modpep_from_sequence(self, peptide_sequence):
        return self._pepname2modpep.get(peptide_sequence)

    def _load_modified_peptides_from_tree(self):
        modified_pepnodes = self._get_modified_peptide_nodes()
        for mod_pep_node in modified_pepnodes:
            self._update_pepname2modpep(mod_pep_node)

    def _get_modified_peptide_nodes(self):
        return anytree.search.findall(self.condpair_tree, lambda x : getattr(x,'type',"") == 'mod_seq', maxlevel=4)

    def _update_pepname2modpep(self, mod_pep_node):
        modified_peptide = PeptideWithSpecificModification(mod_pep_node, self.specific_modification)
        if modified_peptide.specific_modification_found:
            self._pepname2modpep[modified_peptide.peptide_sequence] = modified_peptide


class PeptideWithSpecificModification(OutlierPeptideInfo):
    def __init__(self, node_modpeptide, specific_modification= "[Phospho (STY)]"):
        self.protein_name = self._get_protein_name(node_modpeptide)
        self.modified_sequence = node_modpeptide.name
        self.specific_modification_found = self._check_for_specific_modification(specific_modification)
        if not self.specific_modification_found:
            return
        self.peptide_sequence = self._get_peptide_sequence(node_modpeptide)
        self.fc = node_modpeptide.fc
        self.quality_score = self._get_quality_score(node_modpeptide)

    def _check_for_specific_modification(self, specific_modification):
        return specific_modification in self.modified_sequence

    def _get_peptide_sequence(self, node_modpeptide):
        pepnode = aqclustutils.find_node_parent_at_level(node_modpeptide, level='seq')
        return pepnode.name

    def _get_protein_name(self, node_modpeptide):
        pepnode = aqclustutils.find_node_parent_at_level(node_modpeptide, level='gene')
        return pepnode.name

# Cell
import numpy as np
class ComplementedClusterLoader():
    def __init__(self, outlier_peptide_loader, modified_peptide_loader):
        self._outlier_peptides = outlier_peptide_loader.outlier_peptides
        self._modified_peptide_loader = modified_peptide_loader
        self.complemented_clusters = []
        self._find_complemented_clusters()

    def _find_complemented_clusters(self):
        for outlier_peptide in self._outlier_peptides:
            modified_peptide = self._get_modified_peptide(outlier_peptide)
            if modified_peptide is not None:
                self.complemented_clusters.append(ComplementedCluster(outlier_peptide, modified_peptide))

    def _get_modified_peptide(self, outlier_peptide):
        return self._modified_peptide_loader.get_modpep_from_sequence(outlier_peptide.peptide_sequence)


class ComplementedCluster():
    def __init__(self, outlier_peptide, modified_peptide):
        self.outlier_peptide = outlier_peptide
        self.modified_peptide = modified_peptide
        self._add_normfc_to_modpep()

    def has_opposite_regulation(self):
        return np.sign(self.outlier_peptide.protnormed_fc) == -np.sign(self.modified_peptide.protnormed_fc)

    def get_quality_score(self):
        return max(self.outlier_peptide.quality_score, self.modified_peptide.quality_score)

    def get_outlier_quality_score(self):
        return self.outlier_peptide.quality_score

    def get_modpep_quality_score(self):
        return self.modified_peptide.quality_score

    def get_min_abs_normfc(self):
        return min(abs(self.outlier_peptide.protnormed_fc), abs(self.modified_peptide.protnormed_fc))

    def get_max_abs_normfc(self):
        return max(abs(self.outlier_peptide.protnormed_fc), abs(self.modified_peptide.protnormed_fc))

    def get_outlier_abs_normfc(self):
        return abs(self.outlier_peptide.protnormed_fc)

    def get_ptm_abs_normfc(self):
        return abs(self.modified_peptide.protnormed_fc)

    def get_ptm_abs_fc(self):
        return abs(self.modified_peptide.fc)

    def get_number_mainclust_peptides(self):
        return self.outlier_peptide.num_mainclust_peptides

    def _add_normfc_to_modpep(self):
        self.modified_peptide.protein_fc = self.outlier_peptide.protein_fc
        self.modified_peptide._calc_protnormed_fc()



# Cell
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats


class ComplementedClusterEvaluator():
    def __init__(self, complmented_clusters):
        self._complemented_clusters = complmented_clusters
        self._fcs_outliers = None
        self._fcs_modpeps = None
        self._assign_fold_change_lists()


    def compare_regulation_directions(self, ax):
        opposite_regulation_overview = [int(x.has_opposite_regulation()) for x in self._complemented_clusters]
        self._plot_regulation_direction_histogram(ax, opposite_regulation_overview)

    def scatter_fold_changes(self,ax):

        num_opposite = sum([np.sign(x[0])==-np.sign(x[1]) for x in zip(self._fcs_outliers, self._fcs_modpeps)])
        num_same = sum([np.sign(x[0])==np.sign(x[1]) for x in zip(self._fcs_outliers, self._fcs_modpeps)])
        LOGGER.info(f"{num_same} same, {num_opposite} opposite")
        sns.scatterplot(x =self._fcs_outliers, y=self._fcs_modpeps,ax=ax)
        ax.set_xlabel("outliers")
        ax.set_ylabel("modified_peptides")
        self._set_axis_limits(ax)
        self._draw_horizontal_vertical_line(ax)

    def calculate_correlation(self):
        r, p = scipy.stats.pearsonr(self._fcs_outliers, self._fcs_modpeps)
        LOGGER.info(f"pval is {p}")
        return r


    def _assign_fold_change_lists(self):
        self._fcs_outliers = list([x.outlier_peptide.protnormed_fc for x in self._complemented_clusters])
        self._fcs_modpeps = list([x.modified_peptide.protnormed_fc for x in self._complemented_clusters])

    def _set_axis_limits(self,ax):
        all_lims = ax.get_xlim() + ax.get_ylim() #returns tuples with the lims
        most_extreme_val =  max((abs(x) for x in all_lims))
        ax.set_xlim(-most_extreme_val, most_extreme_val)
        ax.set_ylim(-most_extreme_val, most_extreme_val)

    def _draw_horizontal_vertical_line(self, ax):
        ax.hlines(y=0, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], colors='black')
        ax.vlines(x=0, ymin = ax.get_ylim()[0], ymax = ax.get_ylim()[1], colors='black')




    @staticmethod
    def _plot_regulation_direction_histogram(ax, opposite_regulation_overview):
        ax.hist(opposite_regulation_overview)




# Cell
class QuantileFilterer():
    def __init__(self, objects_to_filter, filterconfigs):
        """Filters objects (e.g. ClusterDiffInfos) based on the
        FilterConfigs defined below. The objects need to have
        a function with the identical name as in FilterConfig.property_name
        which returns a scalar. The quantiles are defined with respect to this scalar.

        Args:
            objects_to_filter (_type_): _description_
            filterconfigs (_type_): _description_
        """

        self._objects_to_filter = objects_to_filter
        self._filterconfigs = filterconfigs

    def get_filtered_list_of_objects(self):
        if self._filterconfigs is None:
            return self._objects_to_filter
        else:
            return self._filter_objects()

    def _filter_objects(self):
        individually_filtered = []
        for filterconf in self._filterconfigs:
            filtered_cclusts = set(self._filter_to_property_quantile(filterconf))
            individually_filtered.append(filtered_cclusts)
        return list(set.intersection(*individually_filtered)) #we filter the quantile of the COMPLETE set for
        #every condition and intersect in the end, alternatively one could successively filter the quantiles


    def _filter_to_property_quantile(self, filterconf):#get the quantiles with the best property scores i.e. quality score
        property_sorted = self._sort_objects_to_filter_by_score(filterconf.property_name)

        if filterconf.quantile_starts_at_lowest:
            number_to_retain = int(filterconf.quantile * len(property_sorted))
            return property_sorted[:number_to_retain]
        else:
            number_to_discard = int((1-filterconf.quantile) * len(property_sorted))
            return property_sorted[number_to_discard:]

    def _sort_objects_to_filter_by_score(self, property_name):
        try:
            return self._property_encodes_instance_variable(property_name)
        except:
            return self._property_encodes_function_name(property_name)

    def _property_encodes_instance_variable(self, property_name):
        return sorted(self._objects_to_filter, key = lambda x : getattr(x, property_name))

    def _property_encodes_function_name(self, property_name):
        return sorted(self._objects_to_filter, key = lambda x : getattr(x, property_name)())



class FilterConfig():
    def __init__(self, property_name, quantile, quantile_starts_at_lowest):
        self.property_name = property_name
        self.quantile = quantile
        self.quantile_starts_at_lowest = quantile_starts_at_lowest



# Cell
class ComplementedClusterFilterer(QuantileFilterer):
    def __init__(self, complemented_clusterloader, clusterfilterconfigs):
        super().__init__(complemented_clusterloader.complemented_clusters, clusterfilterconfigs)

    def get_filtered_complemented_clusters(self):
        return self.get_filtered_list_of_objects()



class ComplementedClusterFilterConfigs():
    def __init__(self, min_abs_normfc_quantile = 1, ptm_abs_normfc_quantile = 1, outlier_abs_normfc_quantile = 1, ptm_absfc_quantile = 1,
    modpep_quality_quantile = 1, outlier_quality_quantile = 1, number_mainclustpeps_quantile = 1):
        self.filterconfigs = []
        self._min_abs_normfc_quantile =min_abs_normfc_quantile
        self._ptm_abs_normfc_quantile=ptm_abs_normfc_quantile
        self._number_mainclustpeps_quantile = number_mainclustpeps_quantile
        self._outlier_abs_normfc_quantile= outlier_abs_normfc_quantile
        self._ptm_absfc_quantile=ptm_absfc_quantile
        self._outlier_quality_quantile=outlier_quality_quantile
        self._modpep_quality_quantile = modpep_quality_quantile
        self._number_mainclustpeps_quantile = number_mainclustpeps_quantile
        self._initialize_filter_configs()

    def _initialize_filter_configs(self):
        self.filterconfigs.append(FilterConfig("get_min_abs_normfc", self._min_abs_normfc_quantile, False))
        self.filterconfigs.append(FilterConfig("get_ptm_abs_normfc", self._ptm_abs_normfc_quantile, False))
        self.filterconfigs.append(FilterConfig("get_ptm_abs_fc", self._ptm_absfc_quantile, False))
        self.filterconfigs.append(FilterConfig("get_outlier_abs_normfc", self._outlier_abs_normfc_quantile, False))
        self.filterconfigs.append(FilterConfig("get_outlier_quality_score", self._outlier_quality_quantile, True))
        self.filterconfigs.append(FilterConfig("get_modpep_quality_score", self._modpep_quality_quantile, True))
        self.filterconfigs.append(FilterConfig("get_number_mainclust_peptides", self._number_mainclustpeps_quantile, False))




# Cell

class DiffClusterFilterer(QuantileFilterer):
    def __init__(self, diffclust_list, diffclustfilterconfigs):
        super().__init__(diffclust_list, diffclustfilterconfigs)

    def get_filtered_diffclust_list(self):
        return self.get_filtered_list_of_objects()


class DiffClusterFilterConfig(FilterConfig):
    def __init__(self, fcdiff_quantile = 1, quality_score_quantile = 1, num_mainclust_peptides_quantile = 1, num_outlierclust_peptides_quantile = 1):
        self.filterconfigs = []
        self._fcdiff_quantile = fcdiff_quantile
        self._quality_score_quantile = quality_score_quantile
        self._num_mainclust_peptides_quantile = num_mainclust_peptides_quantile
        self._num_outlierclust_peptides_quantile = num_outlierclust_peptides_quantile
        self._initialize_filter_configs()

    def _initialize_filter_configs(self):
        self.filterconfigs.append(FilterConfig("fcdiff", self._fcdiff_quantile,False))
        self.filterconfigs.append(FilterConfig("quality_score", self._quality_score_quantile, True))
        self.filterconfigs.append(FilterConfig("get_num_mainclust_peptides", self._num_mainclust_peptides_quantile, False))
        self.filterconfigs.append(FilterConfig("get_num_outlierclust_peptides", self._num_outlierclust_peptides_quantile, False))




# Cell

class OutlierPeptideFilterer(QuantileFilterer):
    def __init__(self, outlier_peptide_list, outlierpeptide_filterconfigs):
        super().__init__(outlier_peptide_list, outlierpeptide_filterconfigs)

    def get_filtered_outlier_peptide_list(self):
        return self.get_filtered_list_of_objects()


class OutlierPeptideFilterConfigs(FilterConfig):
    def __init__(self, quality_score_quantile = 1, num_mainclust_peptides_quantile = 1, protnormed_fc_quantile = 1):
        self.filterconfigs = []
        self._quality_score_quantile = quality_score_quantile
        self._num_mainclust_peptides_quantile = num_mainclust_peptides_quantile
        self._protnormed_fc_quantile = protnormed_fc_quantile
        self._initialize_filter_configs()

    def _initialize_filter_configs(self):
        self.filterconfigs.append(FilterConfig("quality_score", self._quality_score_quantile, True))
        self.filterconfigs.append(FilterConfig("num_mainclust_peptides", self._num_mainclust_peptides_quantile, False))
        self.filterconfigs.append(FilterConfig("protnormed_fc", self._protnormed_fc_quantile, False))
