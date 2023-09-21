
# Cell
import pandas as pd
import anytree
import alphaquant.diffquant.diffutils as aqdiffutils
import alphaquant.cluster.cluster_utils as aqclustutils
import alphaquant.utils.utils as aqutils
import alphaquant.viz.visualizations as aqviz
import alphaquant.config.variables as aqvars


class CondpairQuantificationInfo():

    def init_from_loaded_tables(self, diffresults_df, normed_df, condpair_root_node, samplemap_df):
        self.diffresults_df = diffresults_df
        self.diffresults_df = self.diffresults_df.set_index('protein')
        self.normed_df = normed_df
        self.condpair_root_node = condpair_root_node
        self.samplemap_df = samplemap_df
        self.sample2cond = self.__get_sample2cond(samplemap_df)
        return self

    def init_by_loading_tables(self, cond1, cond2, results_dir, samplemap):
        diffresults_df = aqviz.get_diffresult_dataframe(cond1, cond2, results_folder= results_dir)
        self.diffresults_df = diffresults_df.set_index("protein")
        self.normed_df = aqviz.get_normed_peptides_dataframe(cond1, cond2, results_folder= results_dir)
        self.condpair_root_node = aqutils.read_condpair_tree(cond1, cond2, results_folder= results_dir)
        self.samplemap_df = aqdiffutils.load_samplemap(samplemap)
        self.sample2cond = self.__get_sample2cond(self.samplemap_df)
        return self

    @staticmethod
    def __get_sample2cond(samplemap_df):
        sample2cond = dict(zip(samplemap_df["sample"], samplemap_df["condition"]))
        return sample2cond





import pandas as pd
class ProteinIntensityDataFrameGetter():

    def __init__(self, quantification_info, ion_header = 'quant_id'):
        self._quantification_info = quantification_info
        self._ion_header = ion_header

    def get_melted_df_all(self, protein_id, specified_level):
        protein_node = self._get_protein_node(protein_id)
        melted_df = ProteinIntensityDfFormatter(self._quantification_info,  protein_id,ion_header = self._ion_header).get_melted_protein_ion_intensity_table()
        melted_df = ProteinQuantDfAnnotator(self._quantification_info, protein_node, specified_level).get_annotated_melted_df(melted_df)
        return melted_df

    def get_melted_df_selected_peptides(self, protein_id, selected_peptides, specified_level):
        melted_df = self.get_melted_df_all(protein_id, specified_level)
        melted_df = melted_df[[x in selected_peptides for x in melted_df["specified_level"]]]
        return melted_df
    
    def get_melted_df_clusterdiffinfo(self, protein_id, clusterdiffinfo, specified_level):
        protein_node = self._get_protein_node(protein_id)
        melted_df = self.get_melted_df_all(protein_id, specified_level)
        melted_df =  ProteinQuantDfProteoformSubsetter(melted_df, protein_node, clusterdiffinfo).subset_melted_df_to_clusterdiffinfo()
        return melted_df

    def get_protein_diffresults(self, protein_id):
        return self._quantification_info.diffresults_df.loc[protein_id]
    
    def _get_protein_node(self, protein_id):
        return anytree.findall_by_attr(self._quantification_info.condpair_root_node, protein_id, maxlevel=2)[0]



class ProteinIntensityDfFormatter():
    def __init__(self, quantification_info, protein_id, ion_header = 'quant_id'):
        self._ion_header = ion_header
        self._protein_id = protein_id
        self._quantification_info = quantification_info


    def get_melted_protein_ion_intensity_table(self):
        samples = self._get_samples_of_condpair()
        protein_df = self._subset_dataframe_to_protein()
        return self._melt_protein_dataframe(protein_df, samples)

    def _get_samples_of_condpair(self):
        return set.intersection(set(self._quantification_info.normed_df.columns), set(self._quantification_info.sample2cond.keys()))

    def _subset_dataframe_to_protein(self):
        return self._quantification_info.normed_df.xs(self._protein_id, level = 0)

    def _melt_protein_dataframe(self, protein_df, samples):
        df_melted = pd.melt(protein_df.reset_index(), value_vars= samples, id_vars=[self._ion_header], value_name="intensity", var_name="sample")
        df_melted["condition"] = [self._quantification_info.sample2cond.get(x) for x in df_melted["sample"]]
        return df_melted
    


import pandas as pd
class ProteinQuantDfAnnotator():

    def __init__(self, quantification_info, protein_node, specified_level):
        self._quantification_info = quantification_info
        self._protein_node = protein_node
        self._specified_level = specified_level

        self._ion2is_included = {}
        self._ion2predscore = {}
        self._ion2level = {}
        self._ion2parent = {}
        self._ion2cluster = {}

    
    def get_annotated_melted_df(self, melted_df):
        IonConsistencyTester.ensure_that_diffresult_ions_are_in_tree_ions(melted_df, self._protein_node)
        self._fill_ion_mapping_dicts()
        return self._annotate_properties_to_melted_df(melted_df)


    def _fill_ion_mapping_dicts(self):
        level_nodes = anytree.findall(self._protein_node, filter_= lambda x : (x.type == self._specified_level))
        for level_node in level_nodes:
            for child in level_node.children:
                for leaf in child.leaves:
                    self._ion2is_included[leaf.name] = aqclustutils.check_if_node_is_included(child)
                    self._ion2predscore[leaf.name] = self._get_predscore_if_possible(child)
                    self._ion2level[leaf.name] = child.name
                    self._ion2parent[leaf.name] = level_node.name
                    self._ion2cluster[leaf.name] = child.cluster

    def _annotate_properties_to_melted_df(self, melted_df):
        melted_df["is_included"] = [self._ion2is_included.get(x, np.nan) for x in melted_df[aqvars.QUANT_ID]]
        melted_df["predscore"] = [self._ion2predscore.get(x, np.nan) for x in melted_df[aqvars.QUANT_ID]]
        melted_df["specified_level"] = [self._ion2level.get(x,np.nan) for x in melted_df[aqvars.QUANT_ID]]
        melted_df["parent_level"] = [self._ion2parent.get(x,np.nan) for x in melted_df[aqvars.QUANT_ID]]
        melted_df["cluster"] = [self._ion2cluster.get(x,np.nan) for x in melted_df[aqvars.QUANT_ID]]

        melted_df = melted_df.dropna(subset=["is_included", "predscore", "specified_level", "cluster"])
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
    def __init__(self, melted_df, protein_node,clusterdiffinfo : ClusterDiffInfo):
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
import anytree


import alphaquant.diffquant.diffutils as aqdiffutils
import anytree

class ProteinClusterPlotter():
    def __init__(self, protein_node, condpair, protein_intensity_df_getter : ProteinIntensityDataFrameGetter,level = 'seq', add_stripplot = False):
        self._protein_node = protein_node
        self._level = level
        self._add_stripplot = add_stripplot
        self._axes = None
        self._fig = None
        self._condpair = condpair
        self._melted_df = self._init_melted_df(protein_intensity_df_getter)

    def plot_all_child_elements(self):
        parent2elements = self._get_parent2elements()
        self._prepare_axes(parent2elements)
        self._label_x_and_y()

        for idx, (_, elements) in enumerate(parent2elements.items()):
            melted_df_subset = self._subset_to_elements(self._melted_df, elements)
            colormap = ClusterColorMapper().get_element2color(melted_df_subset)
            fcplotter = IonFoldChangePlotter(melted_df=melted_df_subset, condpair = self._condpair, add_stripplot=self._add_stripplot)
            fcplotter.plot_fcs_with_specified_color_scheme(colormap,self._axes[idx])
            #self._set_title_of_subplot(ax = self._axes[idx], peptide_nodes = cluster_sorted_groups_of_peptide_nodes[idx], first_subplot=idx==0)
        plt.show()

    
    def _init_melted_df(self, protein_intensity_df_getter):
        return protein_intensity_df_getter.get_melted_df_all(self._protein_node.name,self._level)

    def _get_parent2elements(self):
        return self._melted_df.groupby('parent_level')['specified_level'].apply(lambda x: list(x.unique())).to_dict()
    



    @staticmethod
    def _subset_to_elements(df_melted, elements):
        return df_melted[[x in elements for x in df_melted["specified_level"]]]


    def _prepare_axes(self, parent2elements):
        num_independent_plots = len(parent2elements.keys())
        width_list = [len(x) for x in parent2elements.values()] #adjust width of each subplot according to peptide number
        factor = 0.5
        total_number_of_peptides = sum(width_list)
        self._fig, self._axes = plt.subplots(1, num_independent_plots,figsize = (total_number_of_peptides*factor,10),sharey=True, sharex=False, gridspec_kw={'width_ratios' : width_list}, squeeze=False)
        self._axes = self._axes[0] #the squeeze=False option always returns a 2D array, even if there is only one subplot

    def _load_level_nodes(self):
        all_child_nodes = []
        nodes_at_level =  anytree.findall(self._protein_node, filter_= lambda x : (x.type == self._level))
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
    def __init__(self):
        self._colormap = aqviz.AlphaPeptColorMap().colorlist
    
    def get_element2color(self, melted_df):
        unique_clusters = melted_df['cluster'].unique()

        cluster2color = {}
        
        num_colors = len(self._colormap)
        
        for idx, cluster in enumerate(unique_clusters):
            cluster2color[cluster] = self._colormap[idx % num_colors]
        
        element2color = melted_df.set_index('specified_level')['cluster'].map(cluster2color).to_dict()
        
        return element2color




# Cell
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt


class IonFoldChangeCalculator():
    def __init__(self, melted_df, condpair):
        self.melted_df = melted_df
        self._condpair = condpair
        self.__calculate_precursors_and_fcs_from_melted_df()

    def __calculate_precursors_and_fcs_from_melted_df(self):

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
        precfc_tuples = sorted(precfc_tuples, key = lambda x : x[0])
        self.precursors = [x[0] for x in precfc_tuples]
        self.fcs = [x[1] for x in precfc_tuples]


class IonFoldChangePlotter():
    def __init__(self, melted_df, condpair, property_column = "predscore", is_included_column="is_included", add_stripplot = False):

        ionfc_calculated = IonFoldChangeCalculator(melted_df, condpair)
        self._property_column = property_column
        self._is_included_column = is_included_column
        self._add_stripplot = add_stripplot
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
