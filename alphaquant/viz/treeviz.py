

import matplotlib.pyplot as plt
import networkx as nx
import anytree
import re
from matplotlib import gridspec
import matplotlib.pyplot as plt
import alphaquant.cluster.cluster_utils as aqcluster_utils
import alphaquant.viz.visualizations as aqviz
import alphaquant.viz.fcviz as aqfcviz


class TreeSorter():
    def __init__(self, plotconfig, protein):
        self._plotconfig = plotconfig
        self._protein = protein
        self._protein_sequence = None
        self._define_protein_sequence_if_applicable()
        

    def get_sorted_tree(self):
        self._sort_tree(self._protein)
        return self._protein

    def _sort_tree(self, node):
        self._reorder_children(node)
        for child in node.children:
            self._sort_tree(child)

    def _reorder_children(self, parent_node):
        # Get and sort children
        sorted_children = self._get_sorted_children_according_to_plotconfig(parent_node)

        for child in sorted_children:
            child.parent = None  
            child.parent = parent_node

    def _get_sorted_children_according_to_plotconfig(self, protein):
        if self._plotconfig._order_peptides_along_protein_sequence and self._protein.type == 'gene' and self._protein_sequence is not None:
            return aqcluster_utils.get_sorted_peptides_by_position_in_protein_seq(protein, self._protein_sequence)
        
        else:
            return aqcluster_utils.get_sorted_peptides_by_cluster(protein)
        
        # else:
        #     return aqcluster_utils.get_sorted_peptides_by_name(protein)
        
    def _define_protein_sequence_if_applicable(self):
        if self._plotconfig._order_peptides_along_protein_sequence:
            self._protein_sequence = aqcluster_utils.get_protein_sequence(self._protein, self._plotconfig.pyteomics_fasta)



class TreePlotter():
    def __init__(self, protein, parent_level, fig = None, ax = None, plotconfig = aqfcviz.PlotConfig()):
        self.protein = protein

        self.fig = fig
        self.ax = ax

        self._parent_level = parent_level
        self._plotconfig = plotconfig
        
        self._shorten_protein_to_level()
        self._define_fig_and_ax()
        self._create_graph()
    
    def _define_fig_and_ax(self):
        if self.fig is None or self.ax is None:
            axis_creator = TreePlotAxisCreator(self.protein, self._plotconfig)
            axis_creator.define_tree_fig_and_ax()
            self.fig = axis_creator.fig
            self.ax = axis_creator.ax_tree

    def _shorten_protein_to_level(self):
        self.protein = aqcluster_utils.shorten_root_to_level(self.protein, self._parent_level)


    def _create_graph(self):
        GraphCreator(self.protein, self.ax, self._plotconfig)




class GraphCreator():
    
    def __init__(self, protein, ax, plotconfig):
        self.graph = nx.DiGraph()
        self._protein = protein
        self._ax = ax
        self._plotconfig = plotconfig
        self._graph_parameters = GraphParameters()
        self._id2anytree_node = dict()
        self._colorlist_hex = None

        self._add_edges(protein)
        self._define_id2anytree_node()
        self._define_colorlist()
        self._format_graph()
    
    def _add_edges(self, protein):
        for child in protein.children:
            self.graph.add_edge(id(protein), id(child))
            self.graph.nodes[id(protein)]['label'] = protein.name_reduced
            self.graph.nodes[id(child)]['label'] = child.name_reduced
            self._add_edges(child)

        
    def _define_id2anytree_node(self):
        for node in anytree.PreOrderIter(self._protein):
            self._id2anytree_node[id(node)] = node

    def _define_colorlist(self):
        self._colorlist_hex = [aqviz.rgb_to_hex(x) for x in self._plotconfig.colorlist]

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
            self._ax.text(x, y, labelstring, verticalalignment='center', horizontalalignment='center', fontsize=10, family='monospace', 
                          weight = "bold", rotation = self._plotconfig.label_rotation)

        nx.draw_networkx_edges(self.graph, pos, ax=self._ax, **self._graph_parameters.edge_options)
    
    def _determine_cluster_color(self, anynode):
        return self._colorlist_hex[anynode.cluster]


    
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

    def __init__(self, protein, plotconfig):
        self.fig = None
        self.ax_tree = None
        self.axes_fcs = None
        
        self._protein = protein
        self._plotconfig = plotconfig

    def define_combined_tree_fc_fig_and_axes(self):
        parent2leaves = aqcluster_utils.get_parent2leaves_dict(self._protein)
        num_independent_plots = len(parent2leaves.keys())
        width_list = [len(x) for x in parent2leaves.values()]
        
        num_leaves = len(self._protein.leaves)
        max_depth = aqcluster_utils.find_max_depth(self._protein)
        
        fig_width = min(max(8, num_leaves * 1.3), 100) * self._plotconfig.rescale_factor_x
        fig_height = max(8, max_depth * 4) * self._plotconfig.rescale_factor_y
        
        self.fig = plt.figure(figsize=(fig_width, fig_height))
        
        small_width = fig_width * self._plotconfig.narrowing_factor_for_fcplot 
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



