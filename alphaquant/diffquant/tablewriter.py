import pandas as pd
import numpy as np
import alphaquant.ptm.phospho_inference as aq_phospho_inference
import anytree
import statsmodels.stats.multitest as mt
import os
import alphaquant.utils.utils as aqutils


class TableFromNodeCreator():
    def __init__(self, condpair_tree, type = "gene", min_num_peptides = 1):
        self.results_df = None

        self._type = type
        self._min_num_peptides = min_num_peptides
        self._condpair_tree = condpair_tree
        self._list_of_nodes = self._get_list_of_nodes()
        self._condpair_name_table = self._get_condpair_name()

        self._define_results_df()
        self._filter_annotate_results_df()

    def _get_list_of_nodes(self):
        return anytree.findall(self._condpair_tree, filter_ = lambda x : x.type == self._type)

    def _get_condpair_name(self):
        return aqutils.get_condpairname(self._condpair_tree.name)

    def _define_results_df(self):
        list_of_dicts = []
        for node in self._list_of_nodes:
            list_of_dicts.append(self._get_node_dict(node))
        self.results_df = pd.DataFrame(list_of_dicts)
        
    def _get_node_dict(self, node):
        typename_dict = {"gene": "protein", "seq": "sequence", "mod_seq" : "modified_sequence"} #map the short name in the node to a more descriptive name. "gene" to "protein" is a bit confusing, I plan to change everything to "gene" in the future
        type_name  = typename_dict.get(self._type, self._type)
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

        if hasattr(node, "summed_intensity"):
            node_dict["summed_intensity"] = node.summed_intensity

        if self._type == "gene":
            node_dict["num_peptides"] = len(node.children)
        
        return node_dict
    
    def _filter_annotate_results_df(self):
        self.results_df = TableAnnotatorFilterer(self.results_df, self._list_of_nodes, self._min_num_peptides).results_df
        self.results_df = QualityScoreNormalizer(self.results_df, self._list_of_nodes).results_df
    

class TableAnnotatorFilterer():

    def __init__(self, results_df, list_of_nodes, min_num_peptides):

        self.results_df = results_df

        self._example_node = list_of_nodes[0]
        self._min_num_peptides = min_num_peptides

        self._filter_annotate_results_df()
    
    def _filter_annotate_results_df(self):
        if self._example_node.type == "gene":
            self.results_df = self._filter_num_peptides()
        self.results_df = self._scatter_pvals()
        self.results_df = self._add_fdr()
    
    def _filter_num_peptides(self):
        return self.results_df[self.results_df["num_peptides"] >= self._min_num_peptides]

    def _scatter_pvals(self): #add some scatter to the pvalues that are 1.00E-16, which is the lowest possible pvalue. This allows for a better visualization as there are less overlapping points. 
        #Scatter is added by adding a very small random number, therefore minimally reducing significance (i.e. not artificially making significance stronger)
        number_of_cut_pvals = (self.results_df['p_value'] == 1.00E-16).sum()
        random_scatter = np.random.uniform(-14.3, -16, size=number_of_cut_pvals)
        random_scatter = 10**random_scatter

        row_has_cut_pval = self.results_df['p_value'] == 1.00E-16
        self.results_df.loc[row_has_cut_pval, 'p_value'] += random_scatter
        return self.results_df



    def _add_fdr(self):
        pvals = self.results_df["p_value"].tolist()
        fdrs = mt.multipletests(pvals, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
        self.results_df["fdr"] = fdrs
        return self.results_df


class RunConfigTableCreator():
    def __init__(self, runconfig):
        self._runconfig = runconfig

        self.runconfig_df = None

        self._define_results_df()

    def _define_results_df(self):
        method_params = self._get_methods_dict_from_runconfig()
        self.runconfig_df = pd.Series(method_params)

    def _get_methods_dict_from_runconfig(self):
        method_params = {}
        local_vars = self._runconfig.__dict__
        for x in local_vars.keys():
            if local_vars[x] is None:
                continue
            if isinstance(local_vars[x], pd.DataFrame):
                continue

            if (("_df" not in x) and ('condpair' not in x) and ('sys'!=x) and ('runconfig' != x)):
                if ("input_file" in x) or ("results_dir" in x):
                    method_params[x] = os.path.abspath(local_vars[x])
                else:
                    method_params[x] = local_vars[x]
        return method_params



class ProteoFormTableCreator():
    def __init__(self, condpair_tree, organism = None):
        self._condpair_tree = condpair_tree
        self._phospho_scorer = PhosphoScorer(organism)
    
        self.proteoform_df = None

        self._define_proteoform_df()
        self._annotate_proteoform_df()
        

    def _define_proteoform_df(self):
        combined_value_dicts = []
        for protein in self._condpair_tree.children:
            value_dict = ValueDictCreator(protein, self._phospho_scorer).value_dict
            combined_value_dicts.append(value_dict)
        combined_dict = self._merge_list_of_dicts(combined_value_dicts)
        self.proteoform_df = pd.DataFrame(combined_dict)
        self.proteoform_df = QualityScoreNormalizer(self.proteoform_df, self._condpair_tree.children[0]).results_df
    
    @staticmethod
    def _merge_list_of_dicts(dict_list):
        combined_dict = {}
        for d in dict_list:
            for key, value in d.items():
                combined_dict.setdefault(key, []).extend(value)
        return combined_dict

    def _annotate_proteoform_df(self):
        self.proteoform_df = ProteoFormTableAnnotator(self.proteoform_df).proteoform_df


class ValueDictCreator():
    def __init__(self, protein, phospho_scorer):

        self._phospho_scorer = phospho_scorer
        self.value_dict = self._get_value_dict_for_protein(protein)
    
    def _get_value_dict_for_protein(self, protein):
        value_dict = {}
        cluster2peptides = self._get_cluster2peptides(protein)
        for cluster, peptides in cluster2peptides.items():
            value_dict["protein"] = value_dict.get("protein", []) + [protein.name]
            value_dict["proteoform_id"] = value_dict.get("proteoform_id", []) + [f"{protein.name}_{cluster}"]
            value_dict["cluster"] = value_dict.get("cluster", []) + [cluster]
            value_dict["is_reference"] = value_dict.get("is_reference", []) + [cluster==0]
            value_dict["peptides"] = value_dict.get("peptides", []) + [self._get_proetoform_peptides(peptides)]
            value_dict["num_peptides"] = value_dict.get("num_peptides", []) + [len(peptides)]
            value_dict["quality_score"] = value_dict.get("quality_score", []) + [self._get_proteoform_quality_score(peptides)]
            value_dict["log2fc"] = value_dict.get("log2fc", []) + [self._get_proteoform_log2fc(peptides)]
            value_dict["fraction_of_peptides"] =  value_dict.get("fraction_of_peptides", []) + [self._get_fraction_of_peptides(peptides, protein)]
            if self._phospho_scorer.phospho_scoring_available:
                value_dict["likely_phospho"] = value_dict.get("likely_phospho", []) + [self._phospho_scorer.check_if_cluster_likely_phospho(peptides)]
        return value_dict

    @staticmethod
    def _get_cluster2peptides(protein):
        cluster2peptides = {}
        for peptide in protein.children:
            cluster2peptides[peptide.cluster] = cluster2peptides.get(peptide.cluster, []) + [peptide]
        return cluster2peptides
    
    @staticmethod
    def _get_proetoform_peptides(peptides):
        return ";".join([peptide.name for peptide in peptides])
    
    def _get_proteoform_quality_score(self, peptides):
        return max([self._get_peptide_quality_score(peptide) for peptide in peptides])#heuristic: take the highest score of all peptides in the proteoform, thereby increasing the chance for a good score, less strong than summation of scores
    
    @staticmethod
    def _get_peptide_quality_score(peptide):
        if hasattr(peptide, "predscore"):
            return peptide.predscore
        else:
            return peptide.fraction_consistent * len(peptide.leaves)
    
    @staticmethod
    def _get_proteoform_log2fc(peptides):
        return np.mean([peptide.fc for peptide in peptides])
    
    @staticmethod
    def _get_fraction_of_peptides(peptides, protein):
        fraction = len(peptides) / len(protein.children)
        return round(fraction, 2)


class PhosphoScorer():
    def __init__(self, organism):
        self._organism = organism
        self._supported_organisms = ["human"]

        self.phospho_scoring_available = False
        self.phospo_peptide_database = None

        self._check_if_scoring_available()
        self._initialize_phospho_peptide_database()

    def _check_if_scoring_available(self):
        if self._organism in self._supported_organisms:
            self.phospho_scoring_available = True
    
    def _initialize_phospho_peptide_database(self):
        if self.phospho_scoring_available:
            self.phospo_peptide_database = aq_phospho_inference.load_dl_predicted_phosphoprone_sequences(organism=self._organism)
        
    def check_if_cluster_likely_phospho(self, peptides):
        number_of_likely_phospho = len(self.phospo_peptide_database.intersection({x.name for x in peptides}))
        fraction_of_likely_phospho = number_of_likely_phospho / len(peptides)
        if fraction_of_likely_phospho > 0.2:
            return True
        else:
            return False


class ProteoFormTableAnnotator():
    def __init__(self, proteoform_df):
        self.proteoform_df = proteoform_df
        self._annotate_fcdiff_column()
    
    def _annotate_fcdiff_column(self):
        all_rows = []
        self.proteoform_df = self.proteoform_df.sort_values(by=["proteoform_id"])
        for protein, group_df in self.proteoform_df.groupby("protein"):
            first_row = group_df.iloc[0]
            ref_fc = first_row["log2fc"]
            for i, row in group_df.iterrows():
                row["fcdiff"] = abs(row["log2fc"] - ref_fc)
                all_rows.append(row)
        self.proteoform_df = pd.DataFrame(all_rows)


class QualityScoreNormalizer():
    def __init__(self,  results_df, example_node):
        self.results_df = results_df
        self._example_node = example_node

        self._normalize_quality_score()
        self._invert_quality_score_if_ml()

    def _normalize_quality_score(self):
        scores = self.results_df['quality_score'].values

        # Z-Score Normalization
        mean = np.mean(scores)
        std_dev = np.std(scores)
        scores_standardized = (scores - mean) / std_dev

        # Min-Max Scaling
        min_val = np.min(scores_standardized)
        max_val = np.max(scores_standardized)
        scores_min_max_scaled = (scores_standardized - min_val) / (max_val - min_val)

        # Assigning the normalized scores back to the DataFrame
        self.results_df['quality_score'] = scores_min_max_scaled

        return self.results_df

    def _invert_quality_score_if_ml(self):
        if hasattr(self._example_node, "predscore"):
            self.results_df["quality_score"] = 1 - self.results_df["quality_score"]
        return self.results_df
