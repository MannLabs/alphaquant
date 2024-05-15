import pandas as pd
import anytree
import alphaquant.cluster.cluster_utils as aqclustutils
import alphaquant.plotting.base_functions as aq_plot_base
import alphaquant.config.variables as aqvars
import alphamap.organisms_data
import alphaquant.utils.utils as aq_utils
import alphaquant.resources.database_loader as aq_db_loader

import alphaquant.config.config as aqconfig
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

class FoldChangeVisualizer():

    def __init__(self, condition1, condition2, results_directory, samplemap_file,
                                                        order_along_protein_sequence = False, organism = 'Human',colorlist = aq_plot_base.ClusterColorMap().colorlist, tree_level = 'seq',
                                                        protein_identifier = 'gene_symbol', label_rotation = 90, add_stripplot = False,
                                                        narrowing_factor_for_fcplot = 1/14, rescale_factor_x = 1.0, rescale_factor_y = 2,
                                                        
                                                        ):
        """
        Class to visualize the peptide fold changes of a protein (precursor, fragment fcs etc an also be visualized). Can be initialized once and subsequently used to visualize different proteins with the visualize_protein function.

        Args:
            condition1 (str): name of the first experimental condition
            condition2 (str): name of the second experimental condition
            results_directory (str): location of AlphaQuant results directory
            samplemap_file (str): path to AlphaQuant samplemap file
            protein_identifier (str): Identifier for proteins. Can be 'gene_symbol' or 'uniprot_id'.
            label_rotation (int): Rotation of x-axis labels.
            add_stripplot (bool): Add stripplot to boxplot.
            narrowing_factor_for_fcplot (float): Factor to narrow the plot.
            rescale_factor_x (float): Rescale factor for x-axis.
            rescale_factor_y (float): Rescale factor for y-axis.
            tree_level (str): Specify which level of the tree to visualize, options are 'seq', 'mod_seq', 'mod_seq_charge', 'ion_type'.
            colorlist (list): List of colors for plotting.
            protein_identifier (str): Identifier for proteins. Can be 'gene_symbol' or 'uniprot_id'.

        """

        self.plotconfig = PlotConfig(label_rotation = label_rotation, add_stripplot = add_stripplot, narrowing_factor_for_fcplot = narrowing_factor_for_fcplot, rescale_factor_x = rescale_factor_x, rescale_factor_y = rescale_factor_y, colorlist = colorlist, protein_identifier = protein_identifier, tree_level = tree_level, organism = organism, order_peptides_along_protein_sequence=order_along_protein_sequence)

        self.quantification_info = CondpairQuantificationInfo((condition1, condition2), results_directory, samplemap_file)

        #load the trees containing the tree-based quantification info
        self.condpair_tree = aq_utils.read_condpair_tree(condition1, condition2, results_folder=results_directory)


        self.protein2node = {x.name : x for x in self.condpair_tree.children}

    def plot_list_of_proteins(self, list_of_proteins):
        """
        Returns:
            list: list of figure objects for the individual proteins.
        """
        results_figures = []
        for protein_of_interest in list_of_proteins:
            if protein_of_interest not in self.protein2node.keys():
                LOGGER.warn(f"Protein {protein_of_interest} not found in the tree.")
                continue
            protein_fig = self.plot_protein(protein_of_interest)
            results_figures.append(protein_fig)
        
        return results_figures
    
    def plot_protein(self, protein_of_interest, selected_peptides = None):
        """
        Returns:
            figure: figure object for the individual protein.
        """
        protein_node = self.protein2node[protein_of_interest]
        cluster_plotter = ProteinPlot(protein_node, self.quantification_info, self.plotconfig, selected_peptides)
        return cluster_plotter.fig




class PlotConfig():
    def __init__(self, label_rotation = 90, add_stripplot = False, narrowing_factor_for_fcplot = 1/14, rescale_factor_x = 1.0, rescale_factor_y = 2, 
                 colorlist = aq_plot_base.AlphaQuantColorMap().colorlist, protein_identifier = 'gene_symbol', tree_level = 'seq', organism = 'Human', 
                 order_peptides_along_protein_sequence = False):
        """
        Configuration class for plotting.

        Attributes:
            protein_identifier (str): Identifier for proteins. Can be 'gene_symbol' or 'uniprot_id'.
            label_rotation (int): Rotation of x-axis labels.
            add_stripplot (bool): Add stripplot to boxplot.
            narrowing_factor_for_fcplot (float): Factor to narrow the plot.
            rescale_factor_x (float): Rescale factor for x-axis.
            rescale_factor_y (float): Rescale factor for y-axis.
            tree_level (str): Specify which level of the tree to visualize, options are 'seq', 'mod_seq', 'mod_seq_charge', 'ion_type'.
            colorlist (list): List of colors for plotting.
            protein_identifier (str): Identifier for proteins. Can be 'gene_symbol' or 'uniprot_id'.
        """
        self.label_rotation = label_rotation
        self.add_stripplot = add_stripplot
        self.narrowing_factor_for_fcplot = narrowing_factor_for_fcplot
        self.rescale_factor_x = rescale_factor_x
        self.rescale_factor_y = rescale_factor_y
        self.colorlist = colorlist
        self.protein_identifier = protein_identifier #can be 'gene_symbol' or 'uniprot_id'
        self.tree_level = tree_level

        self.parent_level = aqclustutils.LEVELS_UNIQUE[aqclustutils.LEVELS_UNIQUE.index(self.tree_level)+1]
        self.order_peptides_along_protein_sequence = order_peptides_along_protein_sequence
        self._order_by_cluster = not order_peptides_along_protein_sequence
        self._organism = organism
        self.protid2seq = None

        if self.order_peptides_along_protein_sequence:
            self._load_sequences()
    
    def _load_sequences(self):
        organism = self._organism.lower()
        if self.protein_identifier == 'gene_symbol':
            self.protid2seq = aq_db_loader.get_genename2sequence_dict(organism)
        else:
            self.protid2seq = aq_db_loader.get_uniprot2sequence_dict(organism)


def get_pyteomics_fasta(organism = 'Human'):
        return alphamap.organisms_data.import_fasta(organism)

class CondpairQuantificationInfo():
    def __init__(self, condpair, results_dir, samplemap):
        """CondpairQuantificationInfo bundles all static information needed for the foldchangeplots
        """
        self.condpair = condpair
        cond1 = condpair[0]
        cond2 = condpair[1]
        self.normed_intensity_df = aq_plot_base.get_normed_peptides_dataframe(cond1, cond2, results_folder= results_dir)
        self.sample2cond = self._get_sample2cond(samplemap)
        self.relevant_samples = self._get_relevant_samples()
        self.diffresults_df = self._get_diffresults_df(cond1, cond2, results_dir)

    def _get_sample2cond(self, samplemap):
        samplemap_df = pd.read_csv(samplemap, sep = "\t")
        sample2cond = dict(zip(samplemap_df["sample"], samplemap_df["condition"]))
        return sample2cond

    def _get_relevant_samples(self):
        relevant_samples = []

        for sample, cond in self.sample2cond.items():
            if cond in self.condpair:
                relevant_samples.append(sample)
        return relevant_samples

    
    def _get_diffresults_df(self, cond1, cond2, results_dir):
        return aq_plot_base.get_diffresult_dataframe(cond1, cond2, results_folder= results_dir).set_index("protein")






from alphaquant.cluster.outlier_scoring import ClusterDiffInfo

class ProteinQuantDfProteoformSubsetter():
    def __init__(self, melted_df, protein_node, clusterdiffinfo : ClusterDiffInfo):
        self._melted_df = melted_df
        self._protein_node = protein_node
        self._clusterdiffinfo = clusterdiffinfo


    def subset_melted_df_to_clusterdiffinfo(self):

        clusterdiff_protein_node = self._get_clusterdiff_protein_node()

        df_melted_reduced = self._reduce_dataframe_to_clusterdiff_ions(clusterdiff_protein_node)

        return df_melted_reduced

    def _get_clusterdiff_protein_node(self):
        return self._clusterdiffinfo.get_clusterdiff_protnode(self._protein_node)

    def _reduce_dataframe_to_clusterdiff_ions(self, clusterdiff_protein_node):
        ions_used = {x.name  for x in clusterdiff_protein_node.leaves}
        return self._melted_df[[x in ions_used for x in self._melted_df[aqvars.QUANT_ID]]]




# Cell
import alphaquant.diffquant.diffutils as aqdiffutils
import alphaquant.plotting.treeutils as aqtreeutils
import anytree


class ProteinPlot():
    def __init__(self, protein_node, quantification_info: CondpairQuantificationInfo, plotconfig : PlotConfig, selected_peptides = None):

        self.fig = None
        self.axes = None

        self._protein_node = protein_node
        self._quantification_info = quantification_info
        self._plotconfig = plotconfig
        self._selected_peptides = selected_peptides

        self._shorten_protein_node_according_to_plotconfig()
        self._subset_protein_node_to_selected_peptides_if_applicable()
        self._sort_tree_according_to_plotconfig()
        self._plot_fcs()
    
    def _shorten_protein_node_according_to_plotconfig(self):
        self._protein_node = aqclustutils.clone_tree(self._protein_node)
        self._protein_node = aqclustutils.shorten_root_to_level(self._protein_node,parent_level=self._plotconfig.parent_level)

    def _subset_protein_node_to_selected_peptides_if_applicable(self): #would also work for different levels at the current implementation, however makes most sense for the peptide level
        if self._selected_peptides is not None:
            peptide_nodes = anytree.findall(self._protein_node, filter_= lambda x : x.level == "seq")
            peptide_nodes_to_exclude = [x for x in peptide_nodes if x.name not in self._selected_peptides]
            for peptide_node in peptide_nodes_to_exclude:
                peptide_node.parent = None

    
    def _sort_tree_according_to_plotconfig(self):
        self._protein_node = aqtreeutils.TreeSorter(self._plotconfig, self._protein_node).get_sorted_tree()
    
    def _plot_fcs(self):
        pcplotter = ProteinClusterPlotter(self._protein_node, self._quantification_info, self._plotconfig)
        self.fig =  pcplotter._fig
        self.axes = pcplotter._axes
    



class ProteinClusterPlotter():
    def __init__(self, protein_node, quantification_info : CondpairQuantificationInfo, plotconfig : PlotConfig,
                 parent2elements = None, fig = None, axes = None):
        
        self._protein_node = protein_node
        self._plotconfig = plotconfig
        self._quantification_info = quantification_info

        self._parent2elements = parent2elements
        self._fig = fig
        self._axes = axes
        self._melted_df = None
        
        self._init_melted_df()
        self._define_parent2elements()
        self._define_fig_and_axes()
        self._plot_all_child_elements()


    def _init_melted_df(self):
        protein_intensity_df_getter = ProteinIntensityDataFrameGetter(self._protein_node, self._quantification_info)
        self._melted_df = protein_intensity_df_getter.get_melted_df_all(self._plotconfig.parent_level)

    def _define_parent2elements(self):# for example you have precursor as a parent and ms1 and ms2 as the leafs
        if self._parent2elements is None:
            self._parent2elements =  aqclustutils.get_parent2leaves_dict(self._protein_node)
    
        
    def _define_fig_and_axes(self):
        if self._fig is None or self._axes is None:
            self._prepare_axes()


    def _plot_all_child_elements(self):

        for idx, (_, elements) in enumerate(self._parent2elements.items()): #each parent is a separate subplot
            
            melted_df_subset = self._subset_to_elements(self._melted_df, elements) 
            colormap = ClusterColorMapper(self._plotconfig.colorlist).get_element2color(melted_df_subset)
            ProteinPlot = IonFoldChangePlotter(melted_df=melted_df_subset, condpair = self._quantification_info.condpair, plotconfig=self._plotconfig)
            ProteinPlot.plot_fcs_with_specified_color_scheme(colormap,self._axes[idx])
            #self._set_title_of_subplot(ax = self._axes[idx], peptide_nodes = cluster_sorted_groups_of_peptide_nodes[idx], first_subplot=idx==0)
        self._set_yaxes_to_same_scale()
        self._set_title()
        

    @staticmethod
    def _subset_to_elements(df_melted, elements):
        return df_melted.set_index("specified_level").loc[elements].reset_index()


    def _prepare_axes(self):
        num_independent_plots = len(self._parent2elements.keys())
        width_list = [len(x) for x in self._parent2elements.values()] #adjust width of each subplot according to peptide number
        total_number_of_peptides = sum(width_list)
        figsize = (total_number_of_peptides*0.5,10)
        self._fig, self._axes = plt.subplots(1, num_independent_plots,figsize = figsize,sharey=True, sharex=False, gridspec_kw={'width_ratios' : width_list}, squeeze=False)
        self._axes = self._axes[0] #the squeeze=False option always returns a 2D array, even if there is only one subplot


    def _set_yaxes_to_same_scale(self):
        min_ylim = min(ax.get_ylim()[0] for ax in self._axes)
        max_ylim = max(ax.get_ylim()[1] for ax in self._axes)
        
        for ax in self._axes:
            ax.set_ylim(min_ylim, max_ylim)
    


    def _sort_parent2elements(self, parent2elements):
        sorted_parent2elements = {}
        
        for parent_name, elements in parent2elements.items():
            
            parent_node = anytree.search.find(self._protein_node, lambda node: node.name == parent_name)
            ordered_children_names = [child.name for child in parent_node.children]
            sorted_elements = sorted(elements, key=lambda x: ordered_children_names.index(x))
            sorted_parent2elements[parent_name] = sorted_elements

        return sorted_parent2elements

    
    def _load_level_nodes(self):
        all_child_nodes = []
        nodes_at_level =  anytree.findall(self._protein_node, filter_= lambda x : (x.type == self._parent_level))
        for node in nodes_at_level:
            all_child_nodes += node.children
        return all_child_nodes

    @staticmethod
    def _get_peptide_names_to_plot(cluster_sorted_groups_of_peptide_nodes, cluster_idx):
        return [x.name for x in cluster_sorted_groups_of_peptide_nodes[cluster_idx]]

    def _get_color_from_list(self, idx):
        modulo_idx = idx % (len(self._colormap)) #if idx becomes larger than the list length, start at 0 again
        return self._colormap[modulo_idx]

    def _label_x_and_y(self):
        self._fig.supylabel("log2FC")

    def _set_title(self):
        self._fig.suptitle(self._protein_node.name)


    def _set_title_of_subplot(self, ax, peptide_nodes, first_subplot):
        title_text = self._get_subplot_title_text(peptide_nodes, first_subplot)
        ax.set_title(title_text)

    def _get_subplot_title_text(self, peptide_nodes, first_subplot):
        median_fc = np.median([x.fc for x in peptide_nodes])
        min_quality_score = min([self._get_quality_score(x) for x in peptide_nodes])
        fc_string = f"{median_fc:.2}"[:4]
        quality_string = f"{min_quality_score:.2}"[:4]
        if first_subplot:
            return f"fc {fc_string}\nquality {quality_string}"
        else:
            return f"{fc_string}\n{quality_string}"

    def _get_quality_score(self, peptide_node):
        has_predscore = hasattr(peptide_node, 'predscore')
        if has_predscore:
            return abs(peptide_node.predscore)
        else:
            return 1/peptide_node.fraction_consistent



import pandas as pd
class ProteinIntensityDataFrameGetter():

    def __init__(self, protein_node, quantification_info : CondpairQuantificationInfo, ion_header = 'quant_id'):
        self._protein_node = protein_node
        self._quantification_info= quantification_info
        self._ion_header = ion_header

    def get_melted_df_all(self, specified_level):
        melted_df = ProteinIntensityDfFormatter( self._protein_node, self._quantification_info, self._ion_header).get_melted_protein_ion_intensity_table()                                                
        melted_df = ProteinQuantDfAnnotator(self._protein_node, specified_level).get_annotated_melted_df(melted_df)
        return melted_df

    def get_melted_df_selected_peptides(self, protein_id, selected_peptides, specified_level):
        melted_df = self.get_melted_df_all(protein_id, specified_level)
        melted_df = melted_df[[x in selected_peptides for x in melted_df["specified_level"]]]
        return melted_df
    
    def get_melted_df_clusterdiffinfo(self, clusterdiffinfo, specified_level):
        melted_df = self.get_melted_df_all(specified_level)
        melted_df =  ProteinQuantDfProteoformSubsetter(melted_df, self._protein_node, clusterdiffinfo).subset_melted_df_to_clusterdiffinfo()
        return melted_df

    def get_protein_diffresults(self, protein_id):
        return self._quantification_info.diffresults_df.loc[protein_id]
    
    def _get_protein_node(self, protein_id):
        return anytree.findall_by_attr(self._quantification_info.condpair_root_node, protein_id, maxlevel=2)[0]


    
class ProteinIntensityDfFormatter():
    def __init__(self, protein_node, quantification_info, ion_header):
        self._protein_node = protein_node
        self._ion_header = ion_header
        self._normed_intensity_df = quantification_info.normed_intensity_df
        self._relevant_samples = quantification_info.relevant_samples
        self._sample2cond = quantification_info.sample2cond
        
    
    def get_melted_protein_ion_intensity_table(self):
        protein_df = self._subset_dataframe_to_protein()
        return self._melt_protein_dataframe(protein_df)
        

    def _subset_dataframe_to_protein(self):
        return self._normed_intensity_df.xs(self._protein_node.name, level = 0)
    
    
    def _melt_protein_dataframe(self, protein_df):
        df_melted = pd.melt(protein_df.reset_index(), value_vars = self._relevant_samples, id_vars=[self._ion_header], value_name="intensity", var_name="sample")
        df_melted["condition"] = [self._sample2cond.get(x) for x in df_melted["sample"]]
        return df_melted
    




import pandas as pd
import re

class ProteinQuantDfAnnotator():

    def __init__(self, protein_node, specified_level):
        self._protein_node = protein_node
        self._specified_level = specified_level

        self._ion2is_included = {}
        self._ion2predscore = {}
        self._ion2level = {}
        self._ion2parent = {}
        self._ion2cluster = {}

    
    def get_annotated_melted_df(self, melted_df):
        IonConsistencyTester.ensure_that_diffresult_ions_are_in_tree_ions(melted_df, self._protein_node)
        self._add_leafname_column(melted_df)
        self._fill_ion_mapping_dicts()
        
        return self._annotate_properties_to_melted_df(melted_df)
    
    def _add_leafname_column(self, melted_df):#in case the tree has been shortened, the names of the leaves 
        #in the tree are not the same as the ones in the melted df and need to be adapted
        parentlevel2regex = {
            "gene": r"(SEQ_[^_]+_)",
            "seq": r"(SEQ_[^_]+_MOD__[^_]+__)",
            "mod_seq": r"(SEQ_[^_]+_MOD__[^_]+__CHARGE_\d+_)",
            "mod_seq_charge": r"(SEQ_[^_]+_MOD__[^_]+__CHARGE_\d+_(?:FRG|MS1))",
            "ion_type": r"(SEQ_.+)"
        }
        if self._specified_level not in parentlevel2regex.keys():
            melted_df["leafname"] = melted_df[aqvars.QUANT_ID]
        
        else:
            pattern = parentlevel2regex[self._specified_level]
            melted_df["leafname"] = [self._get_new_leafname(pattern, x) for x in melted_df[aqvars.QUANT_ID]]


    def _get_new_leafname(self, pattern, base_ion_name):
        match = re.search(pattern, base_ion_name)
        if match:
            return match.group(1)
        else: 
            raise Exception(f"Could not parse {base_ion_name} at level {self._specified_level}")


    def _fill_ion_mapping_dicts(self):
        level_nodes = anytree.findall(self._protein_node, filter_= lambda x : (x.level == self._specified_level))
        for level_node in level_nodes:
            for child in level_node.children:
                for leaf in child.leaves:
                    self._ion2is_included[leaf.name] = aqclustutils.check_if_node_is_included(child)
                    self._ion2predscore[leaf.name] = self._get_predscore_if_possible(child)
                    self._ion2level[leaf.name] = child.name
                    self._ion2parent[leaf.name] = level_node.name
                    self._ion2cluster[leaf.name] = child.cluster

    def _annotate_properties_to_melted_df(self, melted_df):
        melted_df["is_included"] = [self._ion2is_included.get(x, np.nan) for x in melted_df["leafname"]]
        melted_df["predscore"] = [self._ion2predscore.get(x, np.nan) for x in melted_df["leafname"]]
        melted_df["specified_level"] = [self._ion2level.get(x,np.nan) for x in melted_df["leafname"]]
        melted_df["parent_level"] = [self._ion2parent.get(x,np.nan) for x in melted_df["leafname"]]
        melted_df["cluster"] = [self._ion2cluster.get(x,np.nan) for x in melted_df["leafname"]]

        columns_to_check = ["is_included", "predscore", "specified_level", "cluster"]

        rows_with_na = melted_df[melted_df[columns_to_check].isna().any(axis=1)]

        if not rows_with_na.empty:
            print("Rows with NA values in the specified columns:")
            LOGGER.warn("NA values detected in the specified columns.")
            LOGGER.info(rows_with_na)
            melted_df = melted_df.dropna(subset=columns_to_check)
        
        return melted_df

    @staticmethod
    def _get_predscore_if_possible(node):
        try:
            return node.predscore
        except:
            return 1.0
    

class IonConsistencyTester():
    @staticmethod
    def ensure_that_diffresult_ions_are_in_tree_ions(df_melted, protein_node):
        protnode_ions = [x.name for x in protein_node.leaves]
        ions_in_df = set(df_melted[aqvars.QUANT_ID]) - set(protnode_ions)
        if len(ions_in_df)>0:
            Exception("Clustered ions are not entirely contained in  observed ions!")


class ClusterColorMapper():
    def __init__(self, colorlist = aq_plot_base.AlphaPeptColorMap().colorlist):
        self._colorlist = colorlist
    
    def get_element2color(self, melted_df):
        unique_clusters = melted_df['cluster'].unique()

        cluster2color = {}
        
        num_colors = len(self._colorlist)
        
        for idx, cluster in enumerate(unique_clusters):
            cluster2color[cluster] = self._colorlist[idx % num_colors]
        
        element2color = melted_df.set_index('specified_level')['cluster'].map(cluster2color).to_dict()
        
        return element2color




# Cell
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt


class IonFoldChangePlotter():
    def __init__(self, melted_df, condpair, property_column = "predscore", is_included_column="is_included", plotconfig = PlotConfig()):

        ionfc_calculated = IonFoldChangeCalculator(melted_df, condpair)
        self._property_column = property_column
        self._is_included_column = is_included_column
        self._add_stripplot = plotconfig.add_stripplot
        self.precursors = ionfc_calculated.precursors
        self.fcs = ionfc_calculated.fcs
        self._melted_df = ionfc_calculated.melted_df

    def plot_ion_selection_overview(self):
        fig, axs = plt.subplots(2, 2,figsize = self._get_fig_width())
        colorgetter = IonPlotColorGetter(melted_df = self._melted_df, property_column=self._property_column, ion_name_column="specified_level", is_included_column=self._is_included_column)

        colormap_relative_strength_all = colorgetter.get_predscore_relative_strength_colormap(set_nonmainclust_elems_whiter=False)
        self.plot_fcs_with_specified_color_scheme(colormap_relative_strength_all, axs[0][0])

        colormap_relative_strength_mainclust = colorgetter.get_predscore_relative_strength_colormap(set_nonmainclust_elems_whiter=True)
        self.plot_fcs_with_specified_color_scheme(colormap_relative_strength_mainclust, axs[1][0])

        colormap_quantiles_all = colorgetter.get_predscore_quantile_colormap(set_nonmainclust_elems_whiter=False)
        self.plot_fcs_with_specified_color_scheme(colormap_quantiles_all, axs[0][1])

        colormap_quantiles_mainclust = colorgetter.get_predscore_quantile_colormap(set_nonmainclust_elems_whiter=True)
        self.plot_fcs_with_specified_color_scheme(colormap_quantiles_mainclust, axs[1][1])

        axs[0][0].set_xticks([], [])
        axs[0][1].set_xticks([], [])
        return fig

    def plot_fcs_predscore_relative_strength(self, set_nonmainclust_elems_white = True, ax = None):
        if ax is None:
            ax = plt.subplot()
        colorgetter = IonPlotColorGetter(melted_df = self._melted_df, property_column=self._property_column, ion_name_column="specified_level", is_included_column=self._is_included_column)
        colormap_relative_strength_all = colorgetter.get_predscore_relative_strength_colormap(set_nonmainclust_elems_whiter=set_nonmainclust_elems_white)
        self.plot_fcs_with_specified_color_scheme(colormap_relative_strength_all, ax)
        return ax

    def plot_fcs_predscore_unicolor(self, color, ax = None):
        if ax is None:
            ax = plt.subplot()
        colorgetter = IonPlotColorGetter(melted_df = self._melted_df, property_column=self._property_column, ion_name_column="specified_level", is_included_column=self._is_included_column)
        colormap_single_color = colorgetter.get_single_color_colormap(color)
        self.plot_fcs_with_specified_color_scheme(colormap_single_color, ax)
        return ax

    def plot_fcs_with_specified_color_scheme(self, colormap, ax):
        if type(colormap) == type(dict()):
            colormap = {idx: colormap.get(self.precursors[idx]) for idx in range(len(self.precursors))}
        
        if self._add_stripplot:
            self._plot_fcs_with_swarmplot(colormap, ax)
        else:
            self._plot_fcs_with_boxplot(colormap, ax)

        idxs = list(range(len(self.precursors)))
        ax.set_xticks(idxs, labels = self.precursors, rotation = 'vertical')
    
    def _plot_fcs_with_swarmplot(self, colormap, ax):
        sns.stripplot(data = self.fcs, ax=ax, palette=colormap)
        sns.boxplot(data = self.fcs, ax=ax, 
            boxprops=dict(facecolor="none", edgecolor="black"))

    def _plot_fcs_with_boxplot(self, colormap, ax):
        sns.boxplot(data = self.fcs, ax=ax, palette=colormap)

    def _get_fig_width(self):
        num_ions = len(self.precursors)
        return (int(0.7*num_ions), 10)
    
class IonFoldChangeCalculator():
    def __init__(self, melted_df, condpair):
        self.melted_df = melted_df
        self._condpair = condpair
        self._calculate_precursors_and_fcs_from_melted_df()

    def _calculate_precursors_and_fcs_from_melted_df(self):

        multiindex_df = self.melted_df.set_index(["condition", aqvars.QUANT_ID])


        df_c1 = multiindex_df.loc[self._condpair[0]]
        df_c2 = multiindex_df.loc[self._condpair[1]]

        precursor2fcs = {}

        for ion in df_c1.index.intersection(df_c2.index):
            intens_c1 = df_c1.loc[ion]["intensity"]
            intens_c2 = df_c2.loc[ion]["intensity"]
            fcs = [x -y for x,y in itertools.product(intens_c1, intens_c2)]
            precursor = df_c1.loc[ion]["specified_level"][0]
            precursor2fcs[precursor] = precursor2fcs.get(precursor, []) + fcs

        precfc_tuples = [(x, y) for x,y in precursor2fcs.items()]
        self.precursors = [x[0] for x in precfc_tuples]
        self.fcs = [x[1] for x in precfc_tuples]


class IonPlotColorGetter():

    def __init__(self, melted_df, property_column, ion_name_column, is_included_column):
        self._melted_df = melted_df
        self._property_column = property_column
        self._ion_name_column = ion_name_column
        self._is_included_column = is_included_column

        self._color_palette = aq_plot_base.AlphaPeptColorMap().colormap_discrete
        self._sorted_map_df = self.__init_sorted_mapping_df()

    def get_predscore_relative_strength_colormap(self, set_nonmainclust_elems_whiter = True):
        max_val = list(self._sorted_map_df[self._property_column])[-1]
        relative_proportions = [x/max_val for x in self._sorted_map_df[self._property_column]] #the lower the predscore the lower the proportion (low values in rgb tuple means darker color)
        colors_derived = [(0.8*x, 0.8*x, 0.8*x) for x in relative_proportions] #rgb_base_level = (0.6, 0.6, 0.6)
        ion_names = [x for x in self._sorted_map_df[self._ion_name_column]]
        name2color = dict(zip(ion_names, colors_derived))

        if set_nonmainclust_elems_whiter:
            name2color = self.__make_nonmainclust_elems_whiter(name2color)

        return name2color

    def get_predscore_quantile_colormap(self, set_nonmainclust_elems_whiter = True):
        sorted_scores = self._sorted_map_df[self._property_column]
        idx_fifty_percent = self.__get_percentile_idx(sorted_scores, 0.5)
        idx_seventy_percent = self.__get_percentile_idx(sorted_scores, 0.7)

        name2color_fifty = self.__map_ionname_to_color(self._color_palette(2), idx_start = 0, idx_end = idx_fifty_percent)
        name2color_seventy = self.__map_ionname_to_color(self._color_palette(1), idx_start=idx_fifty_percent, idx_end=idx_seventy_percent)
        name2color_rest = self.__map_ionname_to_color(self._color_palette(0), idx_start=idx_seventy_percent, idx_end=len(sorted_scores))
        name2color_all = self.__merge_dictionaries([name2color_fifty, name2color_seventy, name2color_rest])
        if set_nonmainclust_elems_whiter:
            name2color_all = self.__make_nonmainclust_elems_whiter(name2color_all)
        return name2color_all

    def get_single_color_colormap(self, color):
        name2color_all = self.__map_ionname_to_color(color, idx_start = 0, idx_end = len(self._sorted_map_df.index))
        return name2color_all



    def __make_nonmainclust_elems_whiter(self, name2color):
        name2is_included = self.__init_name2is_included_map()

        for name in name2is_included.keys():
            if not name2is_included.get(name):
                modified_color =self.__make_color_whiter(name2color.get(name), factor=1.0)
                name2color[name] = modified_color

        return name2color

    def __make_color_whiter(self, color_rgba, factor =1):
        modified_rgba = []
        for val in color_rgba:
            new_val = val + (1-val)*(factor)
            modified_rgba.append(new_val)

        return tuple(modified_rgba)


    def __init_sorted_mapping_df(self):
        return self._melted_df[[self._ion_name_column, self._property_column]].drop_duplicates().sort_values(by = self._property_column, ascending = True, ignore_index = True)

    def __init_name2is_included_map(self):
        return dict(zip(self._melted_df[self._ion_name_column], self._melted_df[self._is_included_column]))

    def __get_percentile_idx(self,sorted_scores, percentile):
        return int(np.floor(percentile*len(sorted_scores)))


    def __merge_dictionaries(self,dicts):
        merged_dicts = {}
        for dict in dicts:
            merged_dicts.update(dict)
        return merged_dicts

    def __map_ionname_to_color(self, color, idx_start, idx_end):
        df_subset = self._sorted_map_df.iloc[idx_start:idx_end]
        colorvec = [color for x in range(idx_start, idx_end)]
        name2color = dict(zip(df_subset[self._ion_name_column], colorvec))
        return name2color
