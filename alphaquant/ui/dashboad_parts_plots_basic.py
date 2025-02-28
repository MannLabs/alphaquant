import os
import re
import glob
import param
import panel as pn
import pandas as pd
import matplotlib.pyplot as plt
import itertools
from io import StringIO

import alphaquant.plotting.fcviz as aq_plot_fcviz
import alphaquant.plotting.base_functions as aq_plot_base

# If using Plotly in Panel
pn.extension('plotly')

class PlottingTab(param.Parameterized):
    """
    A param-based Panel class for protein visualization with volcano plot and FoldChangeVisualizer.
    """
    # Public Parameters
    results_dir = param.String(doc="Path to the results folder containing *_VS_*.results.tsv files")
    samplemap_file = param.String(doc="Path to the samplemap file")

    # Dynamic widgets
    condpairname_select = param.ClassSelector(class_=pn.widgets.Select)
    protein_input = param.ClassSelector(class_=pn.widgets.AutocompleteInput)
    tree_level_select = param.ClassSelector(class_=pn.widgets.Select)

    # Internals
    cond_pairs = []
    cond1 = None
    cond2 = None
    result_df = pd.DataFrame()
    fc_visualizer = None

    def __init__(self, state, **params):
        super().__init__(**params)
        self.state = state

        # Create widgets
        self.condpairname_select = pn.widgets.Select(
            name="Select Condition Pair",
            options=["No conditions"],
            width=300
        )
        self.condpairname_select.param.watch(self._on_condpair_selected, 'value')

        # Add show plots button
        self.show_plots_button = pn.widgets.Button(
            name="Show Plots for Selected Condition Pair",
            button_type="primary",
            width=150,
            disabled=True
        )
        self.show_plots_button.on_click(self._on_show_plots_clicked)

        self.protein_input = pn.widgets.AutocompleteInput(
            name='Select Protein',
            placeholder="Search protein, or click a point in the volcano plot above",
            min_characters=0,
            disabled=True,
            width=400,
            margin=(5, 5, 5, 5)
        )
        self.protein_input.param.watch(self._on_protein_selected, 'value')

        self.tree_level_select = pn.widgets.Select(
            name="Tree Level",
            options=['seq', 'base'],
            value='seq',
            width=200
        )
        self.tree_level_select.param.watch(self._on_tree_level_changed, 'value')

        # Input fields for paths
        self.results_dir_input = pn.widgets.TextInput(
            name='Results Directory:',
            value=self.results_dir,
            placeholder='Enter path to results directory',
            width=600
        )
        self.results_dir_input.param.watch(self._update_from_input, 'value')

        # Plot panes
        self.volcano_pane = pn.Column()
        self.protein_plot_pane = pn.Column()

        # Create a container for protein selection controls that will be hidden initially
        self.protein_controls = pn.Column(
            self.protein_input,
            pn.Row(self.tree_level_select),
            visible=False  # Hide initially
        )

        # Create a title for the protein detail section
        self.protein_section_title = pn.pane.Markdown(
            "## Protein Detail Plot\nSelect a protein by clicking on a point in the volcano plot above",
            visible=False  # Hide initially, just like protein_controls
        )

        # Construct layout
        self.main_layout = pn.Column(
            "## Protein Visualization",
            self.results_dir_input,
            self.condpairname_select,
            self.show_plots_button,
            self.volcano_pane,
            pn.layout.Divider(),  # Add a divider here for visual separation
            self.protein_section_title,  # Add the section title
            self.protein_controls,  # Use the container instead of individual widgets
            self.protein_plot_pane,
            sizing_mode='stretch_width'
        )

        # Initialize if directory provided
        if self.results_dir:
            self._extract_condpairs()

        # Watch for state changes
        self.state.param.watch(self._on_state_results_dir_changed, 'results_dir')

    def _update_from_input(self, event):
        """Direct handler for input widget changes"""
        if event.new:
            self.results_dir = event.new
            self.param.trigger('results_dir')  # Trigger the parameter change
            self._extract_condpairs()
            self.samplemap_file = os.path.join(event.new, "samplemap.tsv")

    def panel(self):
        """Return the main panel layout."""
        return self.main_layout

    def on_results_dir_changed(self, new_value):
        """Handle changes to results directory from other components.
        !the method name has to follow the naming pattern on_<param>_changed in order to be recognized by the state manager
        """
        if isinstance(new_value, param.Event):
            value = new_value.new
        elif hasattr(new_value, 'new'):  # Handle Panel event objects
            value = new_value.new
        else:
            value = new_value

        if value is not None:
            value = str(value)
            self.results_dir = value
            if self.results_dir_input.value != value:
                self.results_dir_input.value = value
            self._extract_condpairs()
            self.samplemap_file = os.path.join(value, "samplemap.tsv")

    def _on_state_results_dir_changed(self, event):
        """Handle changes to results directory from state."""
        if event.new and event.new != self.results_dir:
            self.results_dir = event.new
            self.results_dir_input.value = event.new
            self._extract_condpairs()
            self.samplemap_file = os.path.join(self.results_dir_input.value, "samplemap.tsv")

    def _extract_condpairs(self):
        """Look for '*_VS_*.results.tsv' in the results_dir and update the condition pairs."""
        self.cond_pairs = []
        if not self.results_dir or not os.path.isdir(self.results_dir):
            self.condpairname_select.options = ["No conditions"]
            return

        pattern = os.path.join(self.results_dir, "*_VS_*.results.tsv")
        files = glob.glob(pattern)

        for f in files:
            basename = os.path.basename(f)
            match = re.match(r'(.*?)_VS_(.*?)\.results\.tsv$', basename)
            if match:
                cond1, cond2 = match.group(1), match.group(2)
                self.cond_pairs.append((cond1, cond2))

        if self.cond_pairs:
            pairs_str = [f"{c1}_VS_{c2}" for c1, c2 in self.cond_pairs]
            self.condpairname_select.options = pairs_str
            # Select first pair by default if available
            if len(pairs_str) > 0:
                self.condpairname_select.value = pairs_str[0]
        else:
            self.condpairname_select.options = ["No conditions"]

    def _on_condpair_selected(self, event):
        """Called when user selects a condition pair."""
        selected_str = event.new
        if selected_str == "No conditions":
            self.cond1 = None
            self.cond2 = None
            self._clear_plots()
            self.show_plots_button.disabled = True
            return

        if "_VS_" not in selected_str:
            return

        self.cond1, self.cond2 = selected_str.split("_VS_")
        if not self.samplemap_file:
            self.samplemap_file = os.path.join(self.results_dir_input.value, "samplemap.tsv")
        self._update_data_for_condpair()
        self.show_plots_button.disabled = False

    def _update_data_for_condpair(self):
        """Load the results data and initialize FoldChangeVisualizer."""
        # Clear existing data
        self.result_df = pd.DataFrame()

        # Load results
        results_file = os.path.join(self.results_dir_input.value, f"{self.cond1}_VS_{self.cond2}.results.tsv")

        if os.path.exists(results_file):
            try:
                self.result_df = aq_plot_base.get_diffresult_dataframe(
                    self.cond1, self.cond2,
                    results_folder=self.results_dir_input.value
                )

                # Initialize FoldChangeVisualizer
                self._update_fc_visualizer()

                # Update protein selector
                if not self.result_df.empty and 'protein' in self.result_df.columns:
                    prot_list = self.result_df['protein'].dropna().unique().tolist()
                    self.protein_input.options = prot_list
                    self.protein_input.disabled = False
                else:
                    self.protein_input.options = []
                    self.protein_input.disabled = True
            except Exception as e:
                self.result_df = pd.DataFrame()

    def _build_volcano_plot(self):
        """Build and display the volcano plot."""
        self.volcano_pane.clear()
        if not self.result_df.empty:
            try:
                volcano_figure = aq_plot_base.plot_volcano_plotly(self.result_df)
                # Enable clicking in the plot configuration
                volcano_figure.update_layout(
                    clickmode='event+select',
                    width=800,  # Set fixed width for volcano plot
                    height=600,
                    title=f"Volcano Plot: {self.cond1} vs {self.cond2}",
                    xaxis_title={"text": f"Log2(FC)<br>{self.cond1} - {self.cond2}"},
                )
                volcano_pane = pn.pane.Plotly(
                    volcano_figure,
                    config={'responsive': True, 'displayModeBar': True},
                    sizing_mode='fixed'  # Changed to fixed size
                )
                # Connect click event
                volcano_pane.param.watch(self._on_volcano_click, 'click_data')
                self.volcano_pane.append(volcano_pane)
            except Exception as e:
                pass

    def _on_volcano_click(self, event):
        """Handle volcano plot click events."""
        if event.new and 'points' in event.new and event.new['points']:
            protein_name = event.new['points'][0]['text']
            self.protein_input.value = protein_name
            # Directly update the protein plot when clicked
            self._update_protein_plot(protein_name)

    def _on_protein_selected(self, event):
        """Handle protein selection."""
        prot_name = event.new
        if prot_name:
            self._update_protein_plot(prot_name)

    def _update_protein_plot(self, protein_name):
        """Update the protein plot using FoldChangeVisualizer."""
        self.protein_plot_pane.clear()

        if self.fc_visualizer:
            # Update tree level
            self.fc_visualizer.plotconfig.tree_level = self.tree_level_select.value
            self.fc_visualizer.plotconfig.figsize = (10, 6)  # Set smaller figure size

            # Generate plot
            fig = self.fc_visualizer.plot_protein(protein_name)

            # Convert matplotlib figure to panel
            self.protein_plot_pane.append(pn.pane.Matplotlib(fig, tight=True))

    def _clear_plots(self):
        """Clear all plots."""
        self.volcano_pane.clear()
        self.protein_plot_pane.clear()
        # Hide protein selection controls and section title
        self.protein_controls.visible = False
        self.protein_section_title.visible = False

    def _update_fc_visualizer(self):
        """Update FoldChangeVisualizer with current settings."""
        if hasattr(self, 'fc_visualizer') and self.cond1 and self.cond2:
            try:
                self.fc_visualizer = aq_plot_fcviz.FoldChangeVisualizer(
                    condition1=self.cond1,
                    condition2=self.cond2,
                    results_directory=self.results_dir_input.value,
                    samplemap_file=self.samplemap_file,
                    tree_level=self.tree_level_select.value
                )
            except Exception as e:
                self.fc_visualizer = None

    def _on_tree_level_changed(self, event):
        """Handle tree level changes.
        !the method name has to follow the naming pattern on_<param>_changed in order to be recognized by the state manager"""
        self._update_fc_visualizer()
        if self.protein_input.value:
            # Update the plot
            self._update_protein_plot(self.protein_input.value)

    def _on_show_plots_clicked(self, event):
        """Handle show plots button click."""
        if self.cond1 and self.cond2:
            self._build_volcano_plot()
            # Show protein selection controls and section title after plots are built
            self.protein_controls.visible = True
            self.protein_section_title.visible = True
