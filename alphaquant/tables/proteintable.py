import pandas as pd
import numpy as np
import anytree
import statsmodels.stats.multitest as mt
import os
import alphaquant.utils.utils as aqutils
import alphaquant.tables.tableutils as aqtableutils


class TableFromNodeCreator():
    def __init__(self, condpair_tree, node_type = "gene", min_num_peptides = 1, annotation_file = None):
        self.results_df = None

        self._node_type = node_type
        self._min_num_peptides = min_num_peptides
        self._annotation_file = annotation_file
        self._condpair_tree = condpair_tree
        self._list_of_nodes = self._get_list_of_nodes()
        self._condpair_name_table = self._get_condpair_name()


        self._define_results_df()
        self._filter_annotate_results_df()

    def _get_list_of_nodes(self):
        return anytree.findall(self._condpair_tree, filter_ = lambda x : x.type == self._node_type)

    def _get_condpair_name(self):
        return aqutils.get_condpairname(self._condpair_tree.name)

    def _define_results_df(self):
        list_of_dicts = []
        for node in self._list_of_nodes:
            list_of_dicts.append(self._get_node_dict(node))
        self.results_df = pd.DataFrame(list_of_dicts)
        
    def _get_node_dict(self, node):
        typename_dict = {"gene" : "protein", "seq" : "sequence", "mod_seq" : "modified_sequence"} #map the short name in the node to a more descriptive name. "gene" to "protein" is a bit confusing, I plan to change everything to "gene" in the future
        type_name  = typename_dict.get(self._node_type, self._node_type)
        node_dict = {}
        node_dict["condition_pair"] = self._condpair_name_table
        node_dict[type_name] = node.name
        node_dict["p_value"] = node.p_val
        node_dict["log2fc"] = node.fc
        node_dict["number_of_ions"] = len(node.leaves)
        if hasattr(node, "predscore"):
            node_dict["quality_score"] = node.predscore
        else:
            node_dict["quality_score"] = node.fraction_consistent * len(node.leaves)

        if hasattr(node, "total_intensity"):
            node_dict["total_intensity"] = node.total_intensity

        if self._node_type == "gene":
            node_dict["num_peptides"] = len(node.children)
        
        return node_dict
    
    def _filter_annotate_results_df(self):
        self.results_df = TableAnnotatorFilterer(self.results_df, self._list_of_nodes, self._min_num_peptides, self._annotation_file).results_df
        self.results_df = aqtableutils.QualityScoreNormalizer(self.results_df, self._list_of_nodes).results_df
    

class TableAnnotatorFilterer():

    def __init__(self, results_df, list_of_nodes, min_num_peptides, annotation_file):

        self.results_df = results_df

        self._level_type = list_of_nodes[0].type
        self._min_num_peptides = min_num_peptides
        self._annotation_file = annotation_file

        self._filter_annotate_results_df()
    
    def _filter_annotate_results_df(self):
        if self._level_type== "gene":
            self._filter_num_peptides()
            self._add_annotation_columns_if_applicable()
        self._scatter_pvals()
        self._add_fdr()
    
    def _filter_num_peptides(self):
        self.results_df[self.results_df["num_peptides"] >= self._min_num_peptides]

    def _add_annotation_columns_if_applicable(self):
        if self._annotation_file is not None:    
            annotation_df = pd.read_csv(self._annotation_file, sep = "\t")
            annotation_df = annotation_df.drop_duplicates(subset = "protein", keep="first")
            self.results_df = self.results_df.merge(annotation_df, on = "protein", how = "left")

    def _scatter_pvals(self): #add some scatter to the pvalues that are 1.00E-16, which we set as the lowest possible pvalue. This allows for a better visualization as there are less overlapping points. 
        #Scatter is added by adding a very small random number, therefore minimally reducing significance (i.e. not artificially making significance stronger)
        rng = np.random.RandomState(123)
        number_of_cut_pvals = (self.results_df['p_value'] == 1.00E-16).sum()
        random_scatter = rng.uniform(-14.3, -16, size=number_of_cut_pvals)
        random_scatter = 10**random_scatter

        row_has_cut_pval = self.results_df['p_value'] == 1.00E-16
        self.results_df.loc[row_has_cut_pval, 'p_value'] += random_scatter

    def _add_fdr(self):
        pvals = self.results_df["p_value"].tolist()
        fdrs = mt.multipletests(pvals, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
        self.results_df["fdr"] = fdrs

