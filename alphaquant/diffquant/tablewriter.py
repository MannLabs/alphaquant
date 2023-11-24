import pandas as pd
import numpy as np

class ProteoFormTableCreator():
    def __init__(self, condpair_tree):
        self.condpair_tree = condpair_tree
    
        self.proteoform_df = None

        self._define_proteoform_df()
        self._annotate_proteoform_df()

    def _define_proteoform_df(self):
        combined_value_dicts = []
        for protein in self.condpair_tree.children:
            combined_value_dicts.append(self._get_value_dict_for_protein(protein))
        combined_dict = self._merge_list_of_dicts(combined_value_dicts)
        self.proteoform_df = pd.DataFrame(combined_dict)

    def _annotate_proteoform_df(self):
        self.proteoform_df = ProteoFormTableAnnotator(self.proteoform_df).proteoform_df
        
        
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
    
    @staticmethod
    def _merge_list_of_dicts(dict_list):
        combined_dict = {}
        for d in dict_list:
            for key, value in d.items():
                combined_dict.setdefault(key, []).extend(value)
        return combined_dict


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
            

