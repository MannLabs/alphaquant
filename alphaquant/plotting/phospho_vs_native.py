import pandas as pd
import alphaquant.plotting.fcviz as aq_plot_fcviz
import alphaquant.cluster.cluster_utils as aq_clus_utils


class CombinedPhosphoAndNativeDataframeGetter(): #"native" means unenriched proteome data
    def __init__(self, protein_node_native, list_of_protein_nodes_phospho, quantification_info_native : aq_plot_fcviz.CondpairQuantificationInfo, 
                 quantification_info_phospho : aq_plot_fcviz.CondpairQuantificationInfo):
        
        self._protein_node_native = protein_node_native
        self._list_of_protein_nodes_phospho = list_of_protein_nodes_phospho
        self._quantification_info_native = quantification_info_native
        self._quantification_info_phospho = quantification_info_phospho


        self.melted_df_combined = None

        self._prepare_nodes()
        self._init_melted_df_combined()

    def _prepare_nodes(self): #the nodes need to be cut to the level of interest, which in this case are the peptides (children of "gene" level)
        self._protein_node_native = aq_clus_utils.clone_tree(self._protein_node_native)
        self._protein_node_native = aq_clus_utils.shorten_root_to_level(self._protein_node_native, parent_level="gene")
        self._list_of_protein_nodes_phospho = [aq_clus_utils.clone_tree(protein_node) for protein_node in self._list_of_protein_nodes_phospho]
        self._list_of_protein_nodes_phospho = [aq_clus_utils.shorten_root_to_level(protein_node, parent_level="gene") for protein_node in self._list_of_protein_nodes_phospho]


    def _init_melted_df_combined(self):
        melted_df_native = self._init_melted_df_native()
        melted_df_native["data_type"] = "native"
        melted_df_phospho = self._init_melted_df_phospho()
        melted_df_phospho["data_type"] = "phospho"

        self.melted_df_combined = pd.concat([melted_df_native, melted_df_phospho], ignore_index=True)


    def _init_melted_df_native(self):
        protein_intensity_df_getter_native = aq_plot_fcviz.ProteinIntensityDataFrameGetter(self._protein_node_native, self._quantification_info_native)
        return protein_intensity_df_getter_native.get_melted_df_all(specified_level="gene")

    def _init_melted_df_phospho(self):
        melted_dfs = []
        for protein_node_phospho in self._list_of_protein_nodes_phospho:
            protein_intensity_df_getter_phospho = aq_plot_fcviz.ProteinIntensityDataFrameGetter(protein_node_phospho, self._quantification_info_phospho)
            melted_df_phospho = protein_intensity_df_getter_phospho.get_melted_df_all(specified_level="gene")
            melted_dfs.append(melted_df_phospho)
        return pd.concat(melted_dfs, ignore_index=True)


