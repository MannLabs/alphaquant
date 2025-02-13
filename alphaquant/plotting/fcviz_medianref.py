import pandas as pd
import matplotlib.pyplot as plt
import os

import glob
import alphaquant.cluster.cluster_utils as aqclustutils
import alphaquant.utils.utils as aq_utils
import alphaquant.plotting.base_functions as aq_plot_base
import alphaquant.plotting.fcviz as aq_plot_fcviz
import alphaquant.config.config as aqconfig
import alphaquant.config.variables as aq_variables
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

class FoldChangeVisualizerMedianref():

    def __init__(self, results_directory, samplemap_file, order_along_protein_sequence = False, organism = 'Human', colorlist = aq_plot_base.AlphaQuantColorMap().colorlist, tree_level = 'seq', protein_identifier = 'gene_symbol', label_rotation = 90, add_stripplot = False, narrowing_factor_for_fcplot = 1/14, rescale_factor_x = 1.0, rescale_factor_y = 2, figsize = None, showfliers = True, condpairs_to_plot = None):

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


        self.plotconfig = aq_plot_fcviz.PlotConfig(label_rotation = label_rotation, add_stripplot = add_stripplot, narrowing_factor_for_fcplot = narrowing_factor_for_fcplot, rescale_factor_x = rescale_factor_x, rescale_factor_y = rescale_factor_y, colorlist = colorlist, protein_identifier = protein_identifier, tree_level = tree_level, organism = organism, order_peptides_along_protein_sequence=order_along_protein_sequence, figsize=figsize, showfliers=showfliers)

        if condpairs_to_plot is None:
            self.condpairs = MedianRefConditionPairGetter(results_directory).condition_pairs
        else:
            self.condpairs = condpairs_to_plot


        samplemap_file = self._get_samplemap_file_w_median_condition(samplemap_file)
        self.quantification_infos = [aq_plot_fcviz.CondpairQuantificationInfo((condition1, condition2), results_directory, samplemap_file) for condition1, condition2 in self.condpairs]

        #load the trees containing the tree-based quantification info
        self.condpair_trees = [aq_utils.read_condpair_tree(condition1, condition2, results_folder=results_directory) for condition1, condition2 in self.condpairs]

        self.protein2nodes = [{x.name : x for x in condpair_tree.children} for  condpair_tree in self.condpair_trees]
        self.protein2peptides_of_interest = self._get_protein2peptides_of_interest(results_directory)

    def plot_protein_over_conditions(self, protein_of_interest, fig = None, axes = None):
            # Get valid condition pairs where protein exists
            valid_pairs = [(idx, pair) for idx, pair in enumerate(self.condpairs)
                        if protein_of_interest in self.protein2nodes[idx]]

            if len(valid_pairs) == 0:
                return None

            fig, axes = self._init_fig_and_axes_if_none(fig, axes)

            return self._plot_individual_figures(protein_of_interest, valid_pairs, fig, axes)

    def _init_fig_and_axes_if_none(self, fig, axes):
        if fig is None or axes is None:
            return plt.subplots(nrows=len(self.condpairs), ncols=1, figsize=(3.5, len(self.condpairs)*1.5), squeeze=False)
        else:
            return fig, axes

    def _get_samplemap_file_w_median_condition(self, samplemap_file):
        return samplemap_file.replace(".tsv", "_w_median.tsv")

    def _get_protein2peptides_of_interest(self, results_directory):
        proteoform_df = pd.read_csv(f"{results_directory}/medianref_proteoforms.tsv", sep="\t") #the peptides mapping to a proteoform/protein are written like "peptide1;peptide2;..."
        protein2peptide_map_df = proteoform_df[["protein", "peptides"]].assign(peptide=proteoform_df["peptides"].str.split(";")).explode("peptide") #converts the proteoform dataframe to a long format dataframe that has one column for protein and one for single peptide
        protein2peptide_map = protein2peptide_map_df.groupby("protein")["peptide"].agg(list).to_dict()
        return protein2peptide_map

    def _plot_individual_figures(self, protein_of_interest, valid_pairs, fig, axes):
        """Helper method to create individual figures (original behavior)"""
        for idx, condpair in valid_pairs:
            protein_node = self.protein2nodes[idx][protein_of_interest]
            try:
                aq_plot_fcviz.ProteinPlot(protein_node, self.quantification_infos[idx], self.plotconfig, selected_peptides=None, fig=fig, axes=axes[idx])
                axes[idx][0].set_title(f"{protein_of_interest} in {condpair[0]} vs {condpair[1]}")
                LOGGER.info(f"Plotted protein {protein_of_interest} for condition pair {condpair}")
            except Exception as e:
                LOGGER.error(f"Failed to plot protein {protein_of_interest} for condition pair {condpair}: {e}")
        return fig, axes

    def _get_available_proteins(self):
        """Get list of all available proteins across all condition pairs."""
        all_proteins = set()
        for protein_dict in self.protein2nodes:
            all_proteins.update(protein_dict.keys())
        return all_proteins


class MedianRefConditionPairGetter():
    def __init__(self, results_directory):
        self.condition_pairs = self._get_condition_pairs(results_directory)

    def _get_condition_pairs(self, results_directory):
        all_files = glob.glob(f"{results_directory}/*_VS_median_reference.normed.tsv")
        all_conds = [os.path.basename(x).split(aq_variables.CONDITION_PAIR_SEPARATOR)[0] for x in all_files]
        return [[x, "median_reference"] for x in all_conds]
