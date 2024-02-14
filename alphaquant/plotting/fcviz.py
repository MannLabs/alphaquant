import pandas as pd
import anytree
import alphaquant.cluster.cluster_utils as aqclustutils
import alphaquant.plotting.base_functions as aqviz
import alphaquant.config.variables as aqvars
import alphamap.organisms_data





class PlotConfig():
    def __init__(self):
        self.label_rotation = 90
        self.add_stripplot = False
        self.narrowing_factor_for_fcplot = 1/14
        self.rescale_factor_x = 1.0
        self.rescale_factor_y = 2
        self.pyteomics_fasta = None
        self.parent_level = 'gene'
        self.colorlist = aqviz.ClusterColorMap().colorlist

        self._order_peptides_along_protein_sequence = False
        self._order_by_cluster = True
    
    def set_config_to_order_along_protein_sequence(self, organism = 'Human'):
        self.pyteomics_fasta = get_pyteomics_fasta(organism)
        self._order_peptides_along_protein_sequence = True
        self._order_by_cluster = False


def get_pyteomics_fasta(organism = 'Human'):
        return alphamap.organisms_data.import_fasta(organism)

class CondpairQuantificationInfo():
    def __init__(self, condpair, results_dir, samplemap):
        """CondpairQuantificationInfo bundles all static information needed for the foldchangeplots
        """
        self.condpair = condpair
        cond1 = condpair[0]
        cond2 = condpair[1]
        self.normed_intensity_df = aqviz.get_normed_peptides_dataframe(cond1, cond2, results_folder= results_dir)
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
        return aqviz.get_diffresult_dataframe(cond1, cond2, results_folder= results_dir).set_index("protein")




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
            print(rows_with_na)
            raise ValueError("NA values detected in the specified columns.")
        
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


class FCPlotter():
    def __init__(self, protein_node, quantification_info: CondpairQuantificationInfo, plotconfig : PlotConfig):

        self.fig = None
        self.axes = None

        self._protein_node = protein_node
        self._quantification_info = quantification_info
        self._plotconfig = plotconfig
        self._shorten_protein_node_according_to_plotconfig()
        self._sort_tree_according_to_plotconfig()
        self._plot_fcs()
    
    def _shorten_protein_node_according_to_plotconfig(self):
        self._protein_node = aqclustutils.clone_tree(self._protein_node)
        self._protein_node = aqclustutils.shorten_root_to_level(self._protein_node,parent_level=self._plotconfig.parent_level)
    
    def _sort_tree_according_to_plotconfig(self):
        self._protein_node = aqtreeutils.TreeSorter(self._plotconfig, self._protein_node).get_sorted_tree()
    
    def _plot_fcs(self):
        pcplotter = ProteinClusterPlotter(self._protein_node, self._quantification_info, self._plotconfig)
        pcplotter.plot_all_child_elements()
        self.fig =  pcplotter._fig
        self.axes = pcplotter._axes
    



class ProteinClusterPlotter():
    def __init__(self, protein_node, quantification_info : CondpairQuantificationInfo, plotconfig : PlotConfig):
        
        self._protein_node = protein_node
        self._plotconfig = plotconfig
        self._quantification_info = quantification_info

        self._fig = None
        self._axes = None
        self._melted_df = None
        
        self._init_melted_df()

    def plot_all_child_elements(self, parent2elements = None, fig = None, axes = None):
        parent2elements = self._get_parent2elements(parent2elements)
        #self._sort_parent2elements(parent2elements)
        self._define_fig_and_axes(fig, axes, parent2elements)

        for idx, (_, elements) in enumerate(parent2elements.items()):
            
            melted_df_subset = self._subset_to_elements(self._melted_df, elements)
            colormap = ClusterColorMapper(self._plotconfig.colorlist).get_element2color(melted_df_subset)
            fcplotter = IonFoldChangePlotter(melted_df=melted_df_subset, condpair = self._quantification_info.condpair, plotconfig=self._plotconfig)
            fcplotter.plot_fcs_with_specified_color_scheme(colormap,self._axes[idx])
            #self._set_title_of_subplot(ax = self._axes[idx], peptide_nodes = cluster_sorted_groups_of_peptide_nodes[idx], first_subplot=idx==0)
        self._set_yaxes_to_same_scale()
        plt.show()
        

    def _init_melted_df(self):
        protein_intensity_df_getter = ProteinIntensityDataFrameGetter(self._protein_node, self._quantification_info)
        self._melted_df = protein_intensity_df_getter.get_melted_df_all(self._plotconfig.parent_level)


    @staticmethod
    def _subset_to_elements(df_melted, elements):
        return df_melted.set_index("specified_level").loc[elements].reset_index()
    
    def _define_fig_and_axes(self, fig, axes, parent2elements):
        if fig is None or axes is None:
            self._prepare_axes(parent2elements)
        else:
            self._fig = fig
            self._axes = axes    


    def _prepare_axes(self, parent2elements):
        num_independent_plots = len(parent2elements.keys())
        width_list = [len(x) for x in parent2elements.values()] #adjust width of each subplot according to peptide number
        total_number_of_peptides = sum(width_list)
        figsize = (total_number_of_peptides*0.5,10)
        self._fig, self._axes = plt.subplots(1, num_independent_plots,figsize = figsize,sharey=True, sharex=False, gridspec_kw={'width_ratios' : width_list}, squeeze=False)
        self._axes = self._axes[0] #the squeeze=False option always returns a 2D array, even if there is only one subplot


    def _set_yaxes_to_same_scale(self):
        min_ylim = min(ax.get_ylim()[0] for ax in self._axes)
        max_ylim = max(ax.get_ylim()[1] for ax in self._axes)
        
        for ax in self._axes:
            ax.set_ylim(min_ylim, max_ylim)
    
    def _get_parent2elements(self, parent2elements):
        if parent2elements is not None:
            return parent2elements
        else:
            return aqclustutils.get_parent2leaves_dict(self._protein_node)


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

class ClusterColorMapper():
    def __init__(self, colorlist = aqviz.AlphaPeptColorMap().colorlist):
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
        fig, axs = plt.subplots(2, 2,figsize = self.__get_fig_width())
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

    def __get_fig_width(self):
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

        self._color_palette = aqviz.AlphaPeptColorMap().colormap_discrete
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
