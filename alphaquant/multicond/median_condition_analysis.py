
import alphaquant.utils.utils as utils
import pandas as pd


def analyze_and_write_median_condition_results(results_dir):
    """
    Analyzes differential analysis results relative to median reference condition results and writes output files.

    Parameters:
    - results_dir (str): The directory containing the result files.

    Outputs:
    - Writes several TSV files to the results directory:
        - medianref_peptides.tsv,  all peptides and their respective fold changes relative to the median reference condition.
        - medianref_proteoforms.tsv, all proteoform groups and their respective fold changes relative to the median reference condition.
        - medianref_protein_avg.tsv, all protein groups and their average fold changes relative to the median reference condition.
        - medianref_protein_alphaquant.tsv, all protein groups and the fold changes of the reference proteoform in each group relative to the median reference condition.
    """
    condpair2tree = MedianRefResultLoader(results_dir).condpair2tree
    peptide_resolved_proteoform_df = MedianRefConditionCombiner(condpair2tree).peptide_resolved_proteoform_df
    combined_proteoform_df_formatter = CombinedProteoformDfFormatter(peptide_resolved_proteoform_df)
    combined_proteoform_df_formatter.peptide_resolved_proteoform_df.to_csv(f"{results_dir}/medianref_peptides.tsv", sep = "\t", index = None)
    combined_proteoform_df_formatter.proteoform_df.to_csv(f"{results_dir}/medianref_proteoforms.tsv", sep = "\t", index = None)
    combined_proteoform_df_formatter.protein_df_average.to_csv(f"{results_dir}/medianref_protein_avg.tsv", sep = "\t", index = None)
    combined_proteoform_df_formatter.protein_df_pform0.to_csv(f"{results_dir}/medianref_protein_alphaquant.tsv", sep = "\t", index = None)

class MedianRefResultLoader():
    """
    Loads and filters result trees from the specified results directory, keeping only those related to median reference conditions.

    Attributes:
    - condpair2tree (Dict[Any, Any]): A dictionary mapping condition pairs to their corresponding result trees.
    """
    def __init__(self, results_dir):
        """
        Parameters:
        - results_dir (str): The directory containing the results files.
        """
        self._results_dir = results_dir
        self.condpair2tree = {}

        self._define_condpair2tree()
    
    def _define_condpair2tree(self):
        self.condpair2tree = utils.read_all_trees_in_results_folder(self._results_dir)
        self._remove_all_condpairs_w_no_medianref()
    
    def _remove_all_condpairs_w_no_medianref(self):
        keys_to_remove = [x for x in self.condpair2tree.keys() if "median_reference" not in x]
        for key in keys_to_remove:
            del self.condpair2tree[key]
    
    

class MedianRefConditionCombiner():
    """
    Combines result trees from different median reference conditions into a single peptide-resolved proteoform dataframe.

    Attributes:
    - peptide_resolved_proteoform_df (pd.DataFrame): The combined dataframe with peptide-resolved proteoforms.
    """
    def __init__(self, condpair2tree):
        self._condpair2tree = condpair2tree

        self.protein2nodes = {}
        self.peptide_resolved_proteoform_df: pd.DataFrame = None

        self._define_protein2nodes()      
        self._define_combined_dataframe()

    def _define_protein2nodes(self):
        """
        Organizes protein nodes from all conditions into a dictionary grouped by protein names.
        """
        for tree in self._condpair2tree.values():
            cond = tree.name[0]
            for protnode in tree.children:
                self._extend_protein2nodes(protnode, cond)

    def _extend_protein2nodes(self, protnode, cond):
        protnode.cond = cond
        if protnode.name not in self.protein2nodes:
            self.protein2nodes[protnode.name] = []
        self.protein2nodes[protnode.name].append(protnode)

    def _define_combined_dataframe(self):
        list_of_proteoform_dfs = []
        for protein, nodes_same_protein_different_conditions in self.protein2nodes.items():
            proteoform_df = ProteoformConditionAligner(nodes_same_protein_different_conditions).proteoform_df
            proteoform_df = self._add_pvalue_column(proteoform_df, nodes_same_protein_different_conditions)
            list_of_proteoform_dfs.append(proteoform_df)
        self.peptide_resolved_proteoform_df = pd.concat(list_of_proteoform_dfs)
    
    def _add_pvalue_column(self, proteoform_df, nodes_same_protein_different_conditions):
        bf_corrected_pvalue = min([x.p_val for x in nodes_same_protein_different_conditions])*len(nodes_same_protein_different_conditions)
        bf_corrected_pvalue = min(bf_corrected_pvalue, 1)
        proteoform_df["p_value"] = bf_corrected_pvalue
        return proteoform_df
        
    


class ProteoformConditionAligner():
    """
    Aligns proteoforms across different conditions for a single protein.

    Attributes:
    - proteoform_df (pd.DataFrame): The aligned proteoform dataframe for the protein.
    """
    def __init__(self, nodes_same_protein_different_conditions):
        self._nodes_same_protein_different_conditions = nodes_same_protein_different_conditions

        self.proteoform_df = None

        self._define_proteoform_df()

    def _define_proteoform_df(self):
        info_df_creator = ProteoformPeptideDfCreator(self._nodes_same_protein_different_conditions)
        groups_of_peptide_clusters = self._define_groups_of_peptide_clusters(info_df_creator.peptide_cluster_df)
        protein_name = self._nodes_same_protein_different_conditions[0].name
        proteoform_df_creator = ProteoformDfCreator(groups_of_peptide_clusters, info_df_creator.peptide_fc_df, protein_name)
        self.proteoform_df = proteoform_df_creator.proteoform_df

    def _define_groups_of_peptide_clusters(self, peptide_cluster_df):
        peptide_cluster_df = peptide_cluster_df.fillna(0) #assigning peptides that are not present in a condition to the main cluster 0, meaning they will not open up a new proteoform group unless they are in a different cluster than 0 in another condition.
        groups = peptide_cluster_df.groupby(list(peptide_cluster_df.columns)).groups #every peptide that has the identical cluster assignments over conditions is grouped together (e.g. all peptides with cluster idxs [1, 0, 2, 0] will be grouped together) 
        groups_of_peptide_clusters = [list(values) for key, values in groups.items()]
        return sorted(groups_of_peptide_clusters, key=lambda x: len(x), reverse=True)


class ProteoformPeptideDfCreator():
    """
    Reformats the nodes from the same protein into dataframes that can be further processed. The dataframes contain the peptide names and their fold changes or cluster assignments.

    Attributes:
    - peptide_cluster_df (pd.DataFrame): DataFrame with peptide cluster assignments.
    - peptide_fc_df (pd.DataFrame): DataFrame with peptide fold changes.
    """
    def __init__(self, nodes_same_protein_different_conditions):
        self._nodes_same_protein_different_conditions = nodes_same_protein_different_conditions

        self.peptide_cluster_df = None
        self.peptide_fc_df = None

        self._define_peptide_cluster_df()
        self._define_peptide_fc_df()

    
    def _define_peptide_cluster_df(self):
        self.peptide_cluster_df = self._get_peptide_df("cluster")
    
    def _define_peptide_fc_df(self):
        self.peptide_fc_df = self._get_peptide_df("fc")

    
    def _get_peptide_df(self, attribute): #attribute is fc or cluster.
        list_of_pepfc_series = []
        for cond_idx in range(len(self._nodes_same_protein_different_conditions)):
            protein_node = self._nodes_same_protein_different_conditions[cond_idx]
            peptides_for_cond = protein_node.children
            pepfc_series = pd.Series(index=[x.name for x in peptides_for_cond], data=[getattr(x, attribute) for x in peptides_for_cond])
            pepfc_series.name = protein_node.parent.name[0]
            list_of_pepfc_series.append(pepfc_series)
        return pd.DataFrame(list_of_pepfc_series).T #The resulting dataframe will return NaN values for peptides that are not present in a condition, via the 
    

class ProteoformDfCreator():
    """
    Creates a proteoform DataFrame from groups of peptides and their fold changes.

    Attributes:
    - proteoform_df (pd.DataFrame): DataFrame containing proteoform information.
    """
    def __init__(self, groups_of_peptide_clusters, peptide_fc_df, protein_name):
        self._groups_of_peptide_clusters = groups_of_peptide_clusters
        self._peptide_fc_df = peptide_fc_df
        self._protein_name = protein_name

        self.proteoform_df = None

        self._define_proteoform_df()
    
    def _define_proteoform_df(self):
        list_of_proteoform_rows = []
        for idx in range(len(self._groups_of_peptide_clusters)):
            group_of_peptides = self._groups_of_peptide_clusters[idx]
            df_row = self._get_row_of_proteoform_df(group_of_peptides, idx)
            list_of_proteoform_rows.append(df_row)
        self.proteoform_df = pd.concat(list_of_proteoform_rows, )
        
    def _get_row_of_proteoform_df(self, group_of_peptides, idx):
        row = self._peptide_fc_df.loc[group_of_peptides, :]
        row.insert(0,"proteoform_id",f"{self._protein_name}_{idx}")
        row.insert(0, "protein", self._protein_name)
        row.insert(0, "peptides", ";".join(group_of_peptides))
        return row

        



import numpy as np

class CombinedProteoformDfFormatter():
    """takes the peptide resolved proteoform df an formats it to proteoform and protein level"""
    def __init__(self, peptide_resolved_proteoform_df):
        self.peptide_resolved_proteoform_df = peptide_resolved_proteoform_df
        self.proteoform_df = None
        self.protein_df_average = None
        self.protein_df_pform0 = None

        self._define_protein_df_average()
        self._define_protein_df_pform0()
        self._define_proteoform_df()
        self._drop_pval_from_peptides()
    

    
    def _define_protein_df_average(self):
        self.protein_df_average = self.peptide_resolved_proteoform_df.groupby('protein').mean(numeric_only = True).reset_index()
    
    def _define_protein_df_pform0(self):
        is_first_proteoform = [x.endswith("_0") for x in self.peptide_resolved_proteoform_df["proteoform_id"]]
        self.protein_df_pform0 = self.peptide_resolved_proteoform_df[is_first_proteoform].groupby('protein').mean(numeric_only = True).reset_index()
    
    def _define_proteoform_df(self):
        aggregation_dict1 = {"protein" : 'first', "peptides" : "first"}
        aggregation_dict2 = {x: "mean" for x in self.peptide_resolved_proteoform_df.columns if x not in ["protein", "proteoform_id", "peptides"]}
        aggregation_dict = {**aggregation_dict1, **aggregation_dict2}
        proteoform_df = self.peptide_resolved_proteoform_df.groupby('proteoform_id').agg(aggregation_dict).reset_index()
        proteoform_df = proteoform_df.drop(columns=["p_value"])
        proteoform_df = self._add_proteoform_info_columns(proteoform_df)
        self.proteoform_df = proteoform_df

    def _add_proteoform_info_columns(self, proteoform_df):
        list_of_sub_dfs = []
        for protein, sub_df in proteoform_df.groupby('protein'):
            sub_df.insert(0, 'number_of_peptides', [len(peptides.split(';')) for peptides in sub_df['peptides']])
            sub_df = sub_df.set_index(['proteoform_id', 'peptides', 'number_of_peptides'])
            # Separating the reference proteoform
            ref_proteoform = sub_df.iloc[0, 1:]

            corr_to_ref = [1]
            is_ref = [True]
            for idx in range(1, len(sub_df.index)):
                comparison_proteoform = sub_df.iloc[idx, 1:]
                corr = self._calc_corr_to_reference(ref_proteoform, comparison_proteoform)
                corr_to_ref.append(corr)
                is_ref.append(False)
            
            sub_df.insert(0, 'corr_to_ref', corr_to_ref)
            sub_df.insert(1, 'is_reference', is_ref)
            list_of_sub_dfs.append(sub_df)
        
        return pd.concat(list_of_sub_dfs).reset_index()

    def _calc_corr_to_reference(self, ref_proteoform, comparison_proteoform):
        arr1 = np.array(ref_proteoform)
        arr2 = np.array(comparison_proteoform)
        mask1 = ~np.isnan(arr1)
        mask2 = ~np.isnan(arr2)
        mask = mask1 & mask2
        corr = np.corrcoef(arr1[mask], arr2[mask])[0,1]
        return corr

    def _drop_pval_from_peptides(self):
        self.peptide_resolved_proteoform_df = self.peptide_resolved_proteoform_df.drop(columns=["p_value"])
