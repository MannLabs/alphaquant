import numpy as np
import pandas as pd
import anytree
import alphaquant.utils.utils as aqutils
import alphaquant.resources.database_loader as aq_db_loader


import alphamap.preprocessing
import alphamap.organisms_data
import alphamap.sequenceplot
import alphamap.uniprot_integration
import alphaquant.plotting.fcviz as aq_plot_fc
import alphaquant.plotting.base_functions as aq_plot_base


class AlphaMapVisualizer:
    def __init__(self, condition1, condition2, results_directory, samplemap_file,
            order_along_protein_sequence = True, organism = 'Human',colorlist = aq_plot_base.ClusterColorMap().colorlist, tree_level = 'seq',
            protein_identifier = 'gene_symbol', label_rotation = 90, add_stripplot = False,
            narrowing_factor_for_fcplot = 1/14, rescale_factor_x = 1.0, rescale_factor_y = 2):
        
        self._fc_visualizer = aq_plot_fc.FoldChangeVisualizer(condition1, condition2, results_directory, samplemap_file,
            order_along_protein_sequence = order_along_protein_sequence, organism = organism, colorlist = colorlist, 
            tree_level = tree_level, protein_identifier = protein_identifier, label_rotation = label_rotation, 
            add_stripplot = add_stripplot, narrowing_factor_for_fcplot = narrowing_factor_for_fcplot, 
            rescale_factor_x = rescale_factor_x, rescale_factor_y = rescale_factor_y)
        
        self._colorlist = self._get_colorlist(self._fc_visualizer)
        self._gene2protein_mapper = Gene2ProteinMapper(self._fc_visualizer.plotconfig._organism)
        self._df_generator = AlphaMapDfGenerator(self._fc_visualizer.condpair_tree, self._gene2protein_mapper, 
                                                 self._fc_visualizer.plotconfig._organism, self._colorlist)
    
    def _get_colorlist(self, fc_visualizer):
        colorlist = fc_visualizer.plotconfig.colorlist
        return aq_plot_base.rgba_to_hex(colorlist)

    def visualize_protein(self, protein):
        """returns 2 plots: 1) fold change plot and 2) alphamap sequence plot"""
        swissprot_id = self._gene2protein_mapper.get_swissprot_id_if_gene(protein)
        fc_plot = self._fc_visualizer.plot_protein(protein)

        alphamap_plot =  alphamap.sequenceplot.plot_peptide_traces(self._df_generator.cluster_dfs,
                    name = self._df_generator.cluster_names,
                    protein = swissprot_id,
                    fasta = self._df_generator.seq_fasta,
                    uniprot=self._df_generator.uniprot,
                    selected_features=['CHAIN','DOMAIN','STRUCTURE', 'MOD_RES', 'TOPO_DOM'],
                    uniprot_feature_dict=alphamap.uniprot_integration.uniprot_feature_dict, 
                    uniprot_color_dict=alphamap.sequenceplot.uniprot_color_dict,
                    trace_colors = self._df_generator.colorlist)
        
        return fc_plot, alphamap_plot
        



class AlphaMapDfGenerator:

    def __init__(self, condpair_node, gene2protein_mapper, organism = 'Human', colorlist = []):
        self._condpair_node = condpair_node
        self._gene2protein_mapper = gene2protein_mapper

        self.colorlist = colorlist
        self.cluster_dfs = []
        self.cluster_names = []

        self.seq_fasta = alphamap.organisms_data.import_fasta(organism)
        self.uniprot = alphamap.organisms_data.import_uniprot_annotation(organism)

        self._define_cluster_dfs()
        self._define_cluster_names()
        

    def _define_cluster_dfs(self):
        df_allclust = self._generate_alphamap_input_df_from_proteome_condpair_node(self._condpair_node)
        unique_clusters = sorted(df_allclust['cluster'].astype('int').unique())
        for cluster in unique_clusters:
            df_cluster = df_allclust[df_allclust['cluster'] == cluster].drop(columns=['cluster'])
            df_cluster_formatted = alphamap.preprocessing.format_input_data(df=df_cluster, fasta = self.seq_fasta, modification_exp = r'\[.*?\]')
            self.cluster_dfs.append(df_cluster_formatted)
        
    def _generate_alphamap_input_df_from_proteome_condpair_node(self, condpair_node):
        rows = []
        for protein in condpair_node.children:
            protein_name = self._gene2protein_mapper.get_swissprot_id_if_gene(protein.name)
            peptides = anytree.findall(protein, filter_=lambda node: node.type == 'seq')
            for peptide in peptides:
                naked_sequence = aqutils.cut_trailing_parts_seqstring(peptide.name) # Replace this with aqutils if needed
                rows.append({"all_protein_ids": protein_name, 
                             "modified_sequence": naked_sequence, 
                             "naked_sequence": naked_sequence, 
                             "cluster": peptide.cluster})
        
        df = pd.DataFrame(rows)
        return df
    
    def _define_cluster_names(self):
        self.cluster_names = [f'cluster {i}' for i in range(len(self.cluster_dfs))]


class Gene2ProteinMapper:
    def __init__(self, organism = 'Human'):
        """
        Often the 'protein' id encodes the gene symbol, so in this case, we need to map the gene symbol to a respective protein id and we use the swissprot database for this.
        """
        self._gene2protein_dict = self._generate_gene2protein_dict(organism)
    
    def get_swissprot_id_if_gene(self, protein):
        if protein in self._gene2protein_dict:
            return self._gene2protein_dict[protein]
        else:
            return protein
    
    def _generate_gene2protein_dict(self, organism):
        organism_smallcaps = organism.lower()
        return aq_db_loader.get_genename2swissprot_dict(organism_smallcaps)
    
