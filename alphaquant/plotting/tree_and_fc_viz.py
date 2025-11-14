import alphaquant.plotting.fcviz as aqfcviz
import alphaquant.plotting.treeviz as aqtreeviz
import alphaquant.plotting.treeutils as aqtreeutils
import alphaquant.cluster.cluster_utils as aqcluster_utils


class CombinedTreeAndFCPlotter():
    def __init__(self, protein_node, quantification_info: aqfcviz.CondpairQuantificationInfo, plotconfig = aqfcviz.PlotConfig()):

        self._protein_node = protein_node

        self.fig = None
        self.ax_tree = None
        self.axes_fcs = None


        self._protein_node = protein_node
        self._quantification_info = quantification_info
        self._plotconfig = plotconfig
        self._shorten_protein_to_level()
        self._sort_tree_according_to_plotconfig()
        self._maybe_hide_root_in_tree_for_single_branch()
        self._define_fig_and_ax()
        self._plot_tree()
        self._plot_fcs()
        #self._format_fig()

    def _shorten_protein_to_level(self):
        self._protein_node = aqcluster_utils.clone_tree(self._protein_node)
        self._protein_node = aqcluster_utils.shorten_root_to_level(self._protein_node, self._plotconfig.parent_level)

    def _sort_tree_according_to_plotconfig(self):
        self._protein_node = aqtreeutils.TreeSorter(self._plotconfig, self._protein_node).get_sorted_tree()

    def _maybe_hide_root_in_tree_for_single_branch(self):
        try:
            already_set = getattr(self._plotconfig, 'hide_root_in_tree', False)
        except Exception:
            already_set = False

        if not already_set:
            try:
                if hasattr(self._protein_node, 'children') and len(self._protein_node.children) <= 1:
                    self._plotconfig.hide_root_in_tree = True
            except Exception:
                pass

    def _define_fig_and_ax(self):
        axis_creator = aqtreeviz.TreePlotAxisCreator(self._protein_node, self._plotconfig)
        axis_creator.define_combined_tree_fc_fig_and_axes()
        self.fig = axis_creator.fig
        self.ax_tree = axis_creator.ax_tree
        self.axes_fcs = axis_creator.axes_fcs

    def _plot_tree(self):
        """Plot tree using enhanced GraphCreator when annotations are enabled."""
        if self._plotconfig.show_node_annotations:
            aqtreeviz.AnnotatedGraphCreator(self._protein_node, self.ax_tree, self._plotconfig)
        else:
            aqtreeviz.GraphCreator(self._protein_node, self.ax_tree, self._plotconfig)

    def _plot_fcs(self):
        aqfcviz.ProteinClusterPlotter(self._protein_node, self._quantification_info, self._plotconfig,
                                                  fig=self.fig, axes=self.axes_fcs) #updates the axes elements with the fc plot
#        parent2leaves = aqcluster_utils.get_parent2leaves_dict(self._protein_node)

    def _format_fig(self):
        self.fig.tight_layout()
