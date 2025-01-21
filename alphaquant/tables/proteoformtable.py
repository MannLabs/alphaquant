import pandas as pd
import numpy as np
import alphaquant.ptm.phospho_inference as aq_phospho_inference
import statsmodels.stats.multitest as mt
import alphaquant.tables.tableutils as aqtableutils


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
            if protein.missingval == True:
                continue
            value_dict = ValueDictCreator(protein, self._phospho_scorer).value_dict
            combined_value_dicts.append(value_dict)
        combined_dict = self._merge_list_of_dicts(combined_value_dicts)
        self.proteoform_df = pd.DataFrame(combined_dict)
        self.proteoform_df = aqtableutils.QualityScoreNormalizer(self.proteoform_df).results_df
    
    @staticmethod
    def _merge_list_of_dicts(dict_list):
        combined_dict = {}
        for d in dict_list:
            for key, value in d.items():
                combined_dict.setdefault(key, []).extend(value)
        return combined_dict

    def _annotate_proteoform_df(self):
        self.proteoform_df = ProteoFormTableAnnotator(self.proteoform_df).proteoform_df

class ProteoFormTableAnnotator():
    def __init__(self, proteoform_df):
        self.proteoform_df = proteoform_df
        self._annotate_fcdiff_column()
        self._annotate_fdr_column()
    
    def _annotate_fcdiff_column(self):
        all_rows = []
        self.proteoform_df = self.proteoform_df.sort_values(by=["proteoform_id"])
        for protein, group_df in self.proteoform_df.groupby("protein"):
            first_row = group_df.iloc[0]
            ref_fc = first_row["log2fc"]
            for i, row in group_df.iterrows():
                row["fcdiff"] = row["log2fc"] - ref_fc
                row["abs_fcdiff"] = abs(row["fcdiff"])
                all_rows.append(row)
        self.proteoform_df = pd.DataFrame(all_rows)
    
    def _annotate_fdr_column(self):
        mask_of_outlier_pforms = self.proteoform_df["proteoform_pval"].notna()
        pvals = self.proteoform_df.loc[mask_of_outlier_pforms, "proteoform_pval"].tolist()
        if len(pvals)>0:
            fdrs = mt.multipletests(pvals, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
            self.proteoform_df["proteoform_fdr"] = np.nan
            self.proteoform_df.loc[mask_of_outlier_pforms, "proteoform_fdr"] = fdrs
        else:
            self.proteoform_df["proteoform_fdr"] = np.nan


class ValueDictCreator():
    def __init__(self, protein, phospho_scorer):

        self._phospho_scorer = phospho_scorer
        self.value_dict = self._get_value_dict_for_protein(protein)
    
    def _get_value_dict_for_protein(self, protein):
        value_dict = {}
        cluster2peptides = self._get_cluster2peptides(protein)
        quality_score_name = "ml_score" if hasattr(protein.children[0], "ml_score") else "consistency_score"
        for cluster, peptides in cluster2peptides.items():
            value_dict["protein"] = value_dict.get("protein", []) + [protein.name]
            value_dict["proteoform_id"] = value_dict.get("proteoform_id", []) + [f"{protein.name}_{cluster}"]
            value_dict["cluster"] = value_dict.get("cluster", []) + [cluster]
            value_dict["is_reference"] = value_dict.get("is_reference", []) + [cluster==0]
            value_dict["peptides"] = value_dict.get("peptides", []) + [self._get_proetoform_peptides(peptides)]
            value_dict["num_peptides"] = value_dict.get("num_peptides", []) + [len(peptides)]
            value_dict[quality_score_name] = value_dict.get(quality_score_name, []) + [self._get_proteoform_quality_score(peptides)]
            value_dict["log2fc"] = value_dict.get("log2fc", []) + [self._get_proteoform_log2fc(peptides)]
            value_dict["proteoform_pval"] = value_dict.get("proteoform_pval", []) + [peptides[0].proteoform_pval]
            value_dict["proteoform_fcfc"] = value_dict.get("proteoform_fcfc", []) + [peptides[0].proteoform_fcfc]
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
        return sum([self._get_peptide_quality_score(peptide) for peptide in peptides])
    
    @staticmethod
    def _get_peptide_quality_score(peptide):
        if hasattr(peptide, "ml_score"):
            return peptide.ml_score
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