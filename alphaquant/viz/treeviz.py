

import matplotlib.pyplot as plt
import networkx as nx
import anytree
import re
from matplotlib import gridspec
import matplotlib.pyplot as plt
import alphaquant.cluster.cluster_utils as aqcluster_utils
import alphaquant.viz.visualizations as aqviz
import alphaquant.viz.fcviz as aqfcviz

class CombinedTreeAndFCPlotter():
    def __init__(self, protein, parent_level, condpair, intensitydfgetter, add_stripplot = False, label_rotation = 0):
        self._protein = protein

        self.fig = None
        self.ax_tree = None
        self.axes_fcs = None


        self._parent_level = parent_level
        self._condpair = condpair
        self._intensitydfgetter = intensitydfgetter
        self._add_stripplot = add_stripplot
        self._label_rotation = label_rotation
        self._shorten_protein_to_level()
        self._define_fig_and_ax()
        self._plot_tree()
        self._plot_fcs()
        #self._format_fig()
    
    def _shorten_protein_to_level(self):
        self._protein = aqcluster_utils.clone_tree(self._protein)
        self._protein = aqcluster_utils.shorten_root_to_level(self._protein, self._parent_level)

    def _define_fig_and_ax(self):
        axis_creator = TreePlotAxisCreator(self._protein)
        axis_creator.define_combined_tree_fc_fig_and_axes()
        self.fig = axis_creator.fig
        self.ax_tree = axis_creator.ax_tree
        self.axes_fcs = axis_creator.axes_fcs

    def _plot_tree(self):
        GraphCreator(self._protein, self.ax_tree, self._label_rotation)

    def _plot_fcs(self):
        pcplotter = aqfcviz.ProteinClusterPlotter(self._protein, self._condpair, self._intensitydfgetter, parent_level = self._parent_level, add_stripplot = self._add_stripplot)
        parent2leaves = aqcluster_utils.get_parent2leaves_dict(self._protein)
        pcplotter.plot_all_child_elements(parent2elements= parent2leaves, fig=self.fig, axes=self.axes_fcs)
    
    def _format_fig(self):
        self.fig.tight_layout()


class TreePlotter():
    def __init__(self, protein, parent_level, fig = None, ax = None, label_rotation = 0):
        self.protein = protein

        self.fig = fig
        self.ax = ax

        self._parent_level = parent_level
        self._label_rotation = label_rotation
        
        self._shorten_protein_to_level()
        self._define_fig_and_ax()
        self._create_graph()
    
    def _define_fig_and_ax(self):
        if self.fig is None or self.ax is None:
            axis_creator = TreePlotAxisCreator(self.protein)
            axis_creator.define_tree_fig_and_ax()
            self.fig = axis_creator.fig
            self.ax = axis_creator.ax_tree

    def _shorten_protein_to_level(self):
        self.protein = aqcluster_utils.shorten_root_to_level(self.protein, self._parent_level)


    def _create_graph(self):
        GraphCreator(self.protein, self.ax, self._label_rotation)


class GraphCreator():
    
    def __init__(self, protein, ax, label_rotation = 0):
        self.graph = nx.DiGraph()
        self._protein = protein
        self._ax = ax
        self._label_rotation = label_rotation
        self._graph_parameters = GraphParameters()
        self._id2anytree_node = dict()
        self._colorlist = aqviz.AlphaPeptColorMap().colorlist_hex

        self._add_edges(protein)
        self._define_id2anytree_node()
        self._format_graph()
    
    def _add_edges(self, protein):
        children_sorted = sorted(protein.children, key=lambda x: x.name_reduced)
        for child in children_sorted:
            self.graph.add_edge(id(protein), id(child))
            self.graph.nodes[id(protein)]['label'] = protein.name_reduced
            self.graph.nodes[id(child)]['label'] = child.name_reduced
            self._add_edges(child)
        
    def _define_id2anytree_node(self):
        for node in anytree.PreOrderIter(self._protein):
            self._id2anytree_node[id(node)] = node

    def _format_graph(self):
        pos = nx.drawing.nx_agraph.graphviz_layout(self.graph, **self._graph_parameters.layout_params)

        for node in self.graph.nodes():
            matching_anynode  = self._id2anytree_node[node]
            is_included = matching_anynode.is_included
            if not is_included:
                self._graph_parameters.node_options["alpha"] = self._graph_parameters.alpha_excluded
            self._graph_parameters.node_options["node_color"] = self._determine_cluster_color(matching_anynode)
            nx.draw_networkx_nodes(self.graph, pos, nodelist=[node], ax=self._ax, **self._graph_parameters.node_options)

        label_dict = nx.get_node_attributes(self.graph, 'label')

        for node, (x, y) in pos.items():
            labelstring = label_dict[node]
            labelstring = TreeLabelFormatter.format_label_string(labelstring)
            self._ax.text(x, y, labelstring, verticalalignment='center', horizontalalignment='center', fontsize=10, family='monospace', weight = "bold", rotation = self._label_rotation)

        nx.draw_networkx_edges(self.graph, pos, ax=self._ax, **self._graph_parameters.edge_options)
    
    def _determine_cluster_color(self, anynode):
        return self._colorlist[anynode.cluster]


    
    @staticmethod
    def render_tree(root):
        for pre, _, node in anytree.RenderTree(root):
            print("%s%s" % (pre, node.name))


class GraphParameters():
    def __init__(self):
        self.included_color = "skyblue"
        self.excluded_color = "lightgrey"
        self.alpha_included = 1.0
        self.alpha_excluded = 0.2

        self.node_options = {
            "node_color": self.included_color,
            "node_size": 1500,
            "linewidths": 1,
            "alpha": self.alpha_included,  # default alpha
        }

        self.edge_options = {
            "edge_color": "black",
            "arrows": True,
        }

        self.label_options = {
            "font_size": 10,
            "font_color": "darkred",
            "font_weight": "bold",
        }
        self.layout_params = {
        "prog": "dot",
        "args": f"-Gnodesep={4.0/5} -Granksep={4.0/5}"
    }



class TreeLabelFormatter:
    @classmethod
    def format_label_string(cls, labelstring):
        labelstring = cls._cut_leading_type_classifier(labelstring)
        labelstring = cls._remove_leading_trailing_underscores(labelstring)
        labelstring = cls._replace_w_linebreaks(labelstring)
        return labelstring

    @staticmethod
    def _cut_leading_type_classifier(input_string):
        return re.sub(r'^[a-zA-Z0-9]+_', '', input_string)

    @staticmethod
    def _remove_leading_trailing_underscores(input_string):
        return input_string.strip('_')

    @staticmethod
    def _replace_w_linebreaks(input_string):
        result = input_string.replace('_', '\n')
        result = result.replace('[', '\n')
        result = result.replace(']', '\n')
        return result    


class TreePlotAxisCreator():

    def __init__(self, protein):
        self.fig = None
        self.ax_tree = None
        self.axes_fcs = None
        
        self._protein = protein

    def define_combined_tree_fc_fig_and_axes(self):
        parent2leaves = aqcluster_utils.get_parent2leaves_dict(self._protein)
        num_independent_plots = len(parent2leaves.keys())
        width_list = [len(x) for x in parent2leaves.values()]
        
        num_leaves = len(self._protein.leaves)
        max_depth = aqcluster_utils.find_max_depth(self._protein)
        
        fig_width = min(max(8, num_leaves * 1.3), 100)
        fig_height = max(8, max_depth * 4)
        
        self.fig = plt.figure(figsize=(fig_width, fig_height))
        
        small_width = fig_width/14  
        width_ratios = [small_width] + width_list + [small_width]
        
        gs = gridspec.GridSpec(2, num_independent_plots + 2, height_ratios=[1, 1], width_ratios=width_ratios)
        
        self.ax_tree = plt.subplot(gs[0, :])  
        
        ax_small_left = plt.subplot(gs[1, 0])
        ax_small_left.axis('off')
        
        ax_small_right = plt.subplot(gs[1, -1])
        ax_small_right.axis('off')
        
        self.axes_fcs = [plt.subplot(gs[1, i + 1]) for i in range(num_independent_plots)]  # Shifted to exclude the small subplots


    def define_tree_fig_and_ax(self):
        max_depth = aqcluster_utils.find_max_depth(self._protein)
        num_leaves = len(self._protein.leaves)
        fig_width = min(max(8, num_leaves * 1.3),100) 
        fig_height = max(8, max_depth * 2) 

        self.fig, self.ax_tree = plt.subplots(figsize=(fig_width, fig_height))



