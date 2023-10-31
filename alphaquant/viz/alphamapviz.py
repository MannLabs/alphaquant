import numpy as np
import pandas as pd
import anytree
import alphaquant.utils.utils as aqutils


import alphamap.preprocessing
import alphamap.organisms_data
import alphamap.sequenceplot
import alphamap.uniprot_integration



class AlphaMapVisualizer:
    def __init__(self, condpair_node, organism = 'Human'):
        self._df_generator = AlphaMapDfGenerator(condpair_node, organism)
    
    def visualize_protein(self, protein):
        """returns fig object"""
        return alphamap.sequenceplot.plot_peptide_traces(self._df_generator.cluster_dfs,
                    name = self._df_generator.cluster_names,
                    protein = protein,
                    fasta = self._df_generator.human_fasta,
                    uniprot=self._df_generator.human_uniprot,
                    selected_features=['CHAIN','DOMAIN','STRUCTURE', 'MOD_RES', 'TOPO_DOM'],
                    uniprot_feature_dict=alphamap.uniprot_integration.uniprot_feature_dict, 
                    uniprot_color_dict=alphamap.sequenceplot.uniprot_color_dict)




class AlphaMapDfGenerator:

    def __init__(self, condpair_node, organism = 'Human'):
        self._condpair_node = condpair_node
        self.cluster_dfs = []
        self.cluster_names = []

        self.human_fasta = alphamap.organisms_data.import_fasta(organism)
        self.human_uniprot = alphamap.organisms_data.import_uniprot_annotation(organism)

        self._define_cluster_dfs()
        self._define_cluster_names()
        

    def _define_cluster_dfs(self):
        df_allclust = self._generate_alphamap_input_df_from_proteome_condpair_node(self._condpair_node)
        unique_clusters = sorted(df_allclust['cluster'].astype('int').unique())
        for cluster in unique_clusters:
            df_cluster = df_allclust[df_allclust['cluster'] == cluster].drop(columns=['cluster'])
            df_cluster_formatted = alphamap.preprocessing.format_input_data(df=df_cluster, fasta = self.human_fasta, modification_exp = r'\[.*?\]')
            self.cluster_dfs.append(df_cluster_formatted)
        
    def _generate_alphamap_input_df_from_proteome_condpair_node(self, condpair_node):
        rows = []
        for protein in condpair_node.children:
            peptides = anytree.findall(protein, filter_=lambda node: node.type == 'seq')
            for peptide in peptides:
                naked_sequence = aqutils.cut_trailing_parts_seqstring(peptide.name) # Replace this with aqutils if needed
                rows.append({"all_protein_ids": protein.name, 
                             "modified_sequence": naked_sequence, 
                             "naked_sequence": naked_sequence, 
                             "cluster": peptide.cluster})
        
        df = pd.DataFrame(rows)
        return df
    
    def _define_cluster_names(self):
        self.cluster_names = [f'cluster {i}' for i in range(len(self.cluster_dfs))]