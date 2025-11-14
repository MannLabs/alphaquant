

import matplotlib.pyplot as plt
import networkx as nx
import anytree
import re
from matplotlib import gridspec
import matplotlib.pyplot as plt
import alphaquant.cluster.cluster_utils as aqcluster_utils
import alphaquant.plotting.base_functions as aqviz
import alphaquant.plotting.fcviz as aqfcviz
import pandas as pd




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

        root_id = id(self._protein)
        hide_root = getattr(self._plotconfig, 'hide_root_in_tree', False)

        nodes_to_draw = [n for n in self.graph.nodes()] if not hide_root else [n for n in self.graph.nodes() if n != root_id]
        edges_to_draw = list(self.graph.edges()) if not hide_root else [(u, v) for (u, v) in self.graph.edges() if u != root_id and v != root_id]

        for node in nodes_to_draw:
            matching_anynode  = self._id2anytree_node[node]
            is_included = matching_anynode.is_included
            if not is_included:
                self._graph_parameters.node_options["alpha"] = self._graph_parameters.alpha_excluded
            self._graph_parameters.node_options["node_color"] = self._determine_cluster_color(matching_anynode)
            # Allow overriding node size from plotconfig
            if hasattr(self._plotconfig, 'node_size') and self._plotconfig.node_size is not None:
                self._graph_parameters.node_options["node_size"] = self._plotconfig.node_size
            nx.draw_networkx_nodes(self.graph, pos, nodelist=[node], ax=self._ax, **self._graph_parameters.node_options)

        label_dict = nx.get_node_attributes(self.graph, 'label')

        for node in nodes_to_draw:
            x, y = pos[node]
            matching_anynode = self._id2anytree_node[node]
            labelstring = label_dict[node]
            labelstring = TreeLabelFormatter.format_label_string(labelstring)

            # Optionally remove labels for leaf nodes if requested (label will be in x-ticks)
            if getattr(self._plotconfig, 'remove_leaf_labels_in_tree', False) and len(matching_anynode.children) == 0:
                labelstring = ""

            # Use different rotation for leaf nodes (bottom level) vs other nodes
            rotation = self._plotconfig.label_rotation if len(matching_anynode.children) == 0 else 0

            self._ax.text(x, y, labelstring, verticalalignment='center', horizontalalignment='center', fontsize=self._plotconfig.node_fontsize, family='monospace',
                          weight = "bold", rotation = rotation)

        nx.draw_networkx_edges(self.graph, pos, edgelist=edges_to_draw, ax=self._ax, **self._graph_parameters.edge_options)

        # Add vertical padding to avoid cutting labels at top/bottom and hide axis frame
        try:
            ys = [pos[n][1] for n in nodes_to_draw]
            if len(ys) > 0:
                y_min = min(ys)
                y_max = max(ys)
                y_range = (y_max - y_min) if (y_max - y_min) != 0 else 1.0
                pad = 0.15 * y_range
                self._ax.set_ylim(y_min - pad, y_max + pad)
            self._ax.margins(x=0.05, y=0.1)
        except Exception:
            pass
        self._ax.axis('off')

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
        self.alpha_included = 0.6  # More transparent nodes
        self.alpha_excluded = 0.3  # More transparent excluded nodes

        self.node_options = {
            "node_color": self.included_color,
            "node_size": 1500,
            "linewidths": 1,
            "alpha": self.alpha_included,  # default alpha
        }

        self.edge_options = {
            "edge_color": "#CCCCCC",  # Light grey edges
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
        labelstring = cls._remove_default_ion_suffix(labelstring)
        labelstring = cls._replace_w_linebreaks(labelstring)
        return labelstring

    @staticmethod
    def _cut_leading_type_classifier(input_string):
        return re.sub(r'^[a-zA-Z0-9]+_', '', input_string)

    @staticmethod
    def _remove_leading_trailing_underscores(input_string):
        return input_string.strip('_')

    @staticmethod
    def _remove_default_ion_suffix(input_string):
        """Remove the default '_noloss_1' suffix from fragment ion labels."""
        return input_string.replace('_noloss_1', '')

    @staticmethod
    def _replace_w_linebreaks(input_string):
        result = input_string.replace('_', '\n')
        result = result.replace('[', '\n')
        result = result.replace(']', '\n')
        return result


class AnnotatedTreeLabelFormatter(TreeLabelFormatter):
    """Enhanced tree label formatter that can add statistical annotations to node labels."""

    @classmethod
    def format_label_with_annotations(cls, labelstring, node, plotconfig):
        """Format label with optional statistical annotations based on plotconfig.

        Args:
            labelstring (str): Base label string from the node
            node (anytree.Node): The tree node containing attributes
            plotconfig (PlotConfig): Configuration containing annotation settings

        Returns:
            str: Formatted label string with optional annotations
        """
        base_label = cls.format_label_string(labelstring)

        if not plotconfig.show_node_annotations:
            return base_label

        annotations = cls.get_annotation_lines(node, plotconfig)
        if annotations:
            return base_label + "\n" + "\n".join(annotations)
        return base_label

    @classmethod
    def get_annotation_lines(cls, node, plotconfig):
        """Return a list of formatted annotation strings for the given node."""
        annotations = []
        for attr in plotconfig.node_annotation_attributes:
            if hasattr(node, attr):
                value = getattr(node, attr)

                # Skip None or NaN values
                if value is None or (isinstance(value, (int, float)) and pd.isna(value)):
                    continue

                # Use custom format if available, otherwise fallback
                if attr in plotconfig.node_annotation_formats:
                    try:
                        formatted = plotconfig.node_annotation_formats[attr].format(value)
                    except (ValueError, TypeError):
                        formatted = f"{attr}={value}"
                else:
                    if isinstance(value, (int, float)):
                        formatted = f"{attr}={value:.3g}"
                    else:
                        formatted = f"{attr}={value}"

                annotations.append(formatted)
        return annotations


class AnnotatedGraphCreator(GraphCreator):
    """Enhanced GraphCreator that supports configurable node annotations."""

    def __init__(self, protein, ax, plotconfig):
        # Initialize with the same parameters as the parent class
        super().__init__(protein, ax, plotconfig)

    def _format_graph(self):
        """Override _format_graph to use the enhanced label formatter."""
        pos = nx.drawing.nx_agraph.graphviz_layout(self.graph, **self._graph_parameters.layout_params)

        root_id = id(self._protein)
        hide_root = getattr(self._plotconfig, 'hide_root_in_tree', False)

        nodes_to_draw = [n for n in self.graph.nodes()] if not hide_root else [n for n in self.graph.nodes() if n != root_id]
        edges_to_draw = list(self.graph.edges()) if not hide_root else [(u, v) for (u, v) in self.graph.edges() if u != root_id and v != root_id]

        for node in nodes_to_draw:
            matching_anynode = self._id2anytree_node[node]
            is_included = matching_anynode.is_included
            if not is_included:
                self._graph_parameters.node_options["alpha"] = self._graph_parameters.alpha_excluded
            else:
                self._graph_parameters.node_options["alpha"] = self._graph_parameters.alpha_included

            self._graph_parameters.node_options["node_color"] = self._determine_cluster_color(matching_anynode)
            # Allow overriding node size from plotconfig
            if hasattr(self._plotconfig, 'node_size') and self._plotconfig.node_size is not None:
                self._graph_parameters.node_options["node_size"] = self._plotconfig.node_size
            nx.draw_networkx_nodes(self.graph, pos, nodelist=[node], ax=self._ax, **self._graph_parameters.node_options)

        label_dict = nx.get_node_attributes(self.graph, 'label')

        for node in nodes_to_draw:
            x, y = pos[node]
            matching_anynode = self._id2anytree_node[node]
            labelstring = label_dict[node]

            # Use the enhanced formatter with annotations
            labelstring = AnnotatedTreeLabelFormatter.format_label_with_annotations(
                labelstring, matching_anynode, self._plotconfig
            )

            # Optionally remove base label for leaf nodes but keep annotations
            if getattr(self._plotconfig, 'remove_leaf_labels_in_tree', False) and len(matching_anynode.children) == 0:
                if self._plotconfig.show_node_annotations:
                    annotations_only = AnnotatedTreeLabelFormatter.get_annotation_lines(matching_anynode, self._plotconfig)
                    labelstring = "\n".join(annotations_only) if annotations_only else ""
                else:
                    labelstring = ""

            # Adjust font size based on annotation content
            fontsize = self._plotconfig.node_fontsize
            if self._plotconfig.show_node_annotations and len(labelstring.split('\n')) > 2:
                fontsize = max(8, self._plotconfig.node_fontsize - 2)

            # Use different rotation for leaf nodes (bottom level) vs other nodes
            rotation = self._plotconfig.label_rotation if len(matching_anynode.children) == 0 else 0

            self._ax.text(x, y, labelstring, verticalalignment='center', horizontalalignment='center',
                          fontsize=fontsize, family='monospace', weight="bold",
                          rotation=rotation)

        nx.draw_networkx_edges(self.graph, pos, edgelist=edges_to_draw, ax=self._ax, **self._graph_parameters.edge_options)

        # Add vertical padding to avoid cutting labels at top/bottom and hide axis frame
        try:
            ys = [pos[n][1] for n in nodes_to_draw]
            if len(ys) > 0:
                y_min = min(ys)
                y_max = max(ys)
                y_range = (y_max - y_min) if (y_max - y_min) != 0 else 1.0
                pad = 0.15 * y_range
                self._ax.set_ylim(y_min - pad, y_max + pad)
            self._ax.margins(x=0.05, y=0.1)
        except Exception:
            pass
        self._ax.axis('off')


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

        # Determine figure size: honor explicit figsize if provided
        if self._plotconfig.figsize is not None:
            fig_width, fig_height = self._plotconfig.figsize
        else:
            num_leaves = len(self._protein.leaves)
            max_depth = aqcluster_utils.find_max_depth(self._protein)
            fig_width = min(max(8, num_leaves * 1.3), 100) * self._plotconfig.rescale_factor_x
            fig_height = max(8, max_depth * 4) * self._plotconfig.rescale_factor_y

        self.fig = plt.figure(figsize=(fig_width, fig_height))

        small_width = fig_width * self._plotconfig.narrowing_factor_for_fcplot
        width_ratios = [small_width] + width_list + [small_width]

        # Use height ratio from plotconfig (tree : fold_change)
        height_ratios = [self._plotconfig.tree_to_fc_height_ratio, 1]
        gs = gridspec.GridSpec(2, num_independent_plots + 2, height_ratios=height_ratios, width_ratios=width_ratios, hspace=self._plotconfig.subplot_spacing)

        self.ax_tree = plt.subplot(gs[0, :])

        ax_small_left = plt.subplot(gs[1, 0])
        ax_small_left.axis('off')

        ax_small_right = plt.subplot(gs[1, -1])
        ax_small_right.axis('off')

        self.axes_fcs = [plt.subplot(gs[1, i + 1]) for i in range(num_independent_plots)]  # Shifted to exclude the small subplots


    def define_tree_fig_and_ax(self):
        # Determine figure size: honor explicit figsize if provided
        if self._plotconfig.figsize is not None:
            self.fig, self.ax_tree = plt.subplots(figsize=self._plotconfig.figsize)
        else:
            max_depth = aqcluster_utils.find_max_depth(self._protein)
            num_leaves = len(self._protein.leaves)
            fig_width = min(max(8, num_leaves * 1.3),100)
            fig_height = max(8, max_depth * 2)
            self.fig, self.ax_tree = plt.subplots(figsize=(fig_width, fig_height))



