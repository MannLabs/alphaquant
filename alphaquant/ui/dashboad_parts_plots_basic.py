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

        self.protein_input = pn.widgets.AutocompleteInput(
            name='Select Protein',
            placeholder="Type to search protein...",
            min_characters=2,
            disabled=True,
            width=400,
            margin=(5, 5, 5, 5)
        )
        self.protein_input.param.watch(self._on_protein_selected, 'value')

        self.tree_level_select = pn.widgets.Select(
            name="Tree Level",
            options=['seq', 'mod_seq', 'mod_seq_charge', 'ion_type', 'base'],
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
        self.results_dir_input.param.watch(self.on_results_dir_changed, 'value')

        # Add file upload widget
        self.samplemap_fileupload = pn.widgets.FileInput(
            name='Upload Sample Map',
            accept=",".join(".tsv", ".csv", ".txt")
            margin=(5, 5, 10, 20)
        )
        self.samplemap_fileupload.param.watch(self._handle_samplemap_upload, 'value')

        # Create a row for samplemap controls
        self.samplemap_controls = pn.Row(
            pn.Column(
                pn.widgets.StaticText(
                    name='',
                    value='Sample Map: No sample map loaded',
                    styles={'font-weight': 'normal'}
                ),
                self.samplemap_fileupload,
            ),
            margin=(5, 5, 5, 5)
        )

        # Plot panes
        self.volcano_pane = pn.Column()
        self.protein_plot_pane = pn.Column()

        # Construct layout
        self.main_layout = pn.Column(
            "## Protein Visualization",
            self.results_dir_input,
            self.samplemap_controls,
            self.condpairname_select,
            self.volcano_pane,
            pn.Row(self.tree_level_select),
            self.protein_input,
            self.protein_plot_pane,
            sizing_mode='stretch_width'
        )

        # Initialize if directory provided
        self._extract_condpairs()

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
            if self.results_dir_input.value != value:
                self.results_dir_input.value = value
                self._extract_condpairs()

    def _handle_samplemap_upload(self, event):
        """Handle new samplemap file uploads."""
        if not event.new:
            return

        try:
            # Determine file type and separator
            file_ext = os.path.splitext(self.samplemap_fileupload.filename)[-1].lower()
            sep = ',' if file_ext == '.csv' else '\t'

            # Parse the uploaded file into DataFrame
            df = pd.read_csv(
                StringIO(event.new.decode('utf-8')),
                sep=sep,
                dtype=str
            )

            # Update the state with the new DataFrame
            self.state.samplemap_df = df
            self.state.notify_subscribers('samplemap_df')

        except Exception as e:
            self.samplemap_controls[0][0].value = f"Error loading sample map: {str(e)}"

    def on_samplemap_df_changed(self, new_df):
        """Handle changes to samplemap DataFrame from other components.
        !the method name has to follow the naming pattern on_<param>_changed in order to be recognized by the state manager
        """
        if not new_df.empty:
            # Update status
            num_samples = len(new_df)
            num_conditions = len(new_df['condition'].unique()) if 'condition' in new_df.columns else 0
            status_text = f"Sample Map: Already loaded {num_samples} samples, {num_conditions} conditions"
            self.samplemap_controls[0][0].value = status_text

            # Update condition pairs and other visualizations
            self._update_condition_pairs_from_df(new_df)
        else:
            self.samplemap_controls[0][0].value = "Sample Map: No sample map loaded"

    def _update_condition_pairs_from_df(self, df):
        """Update condition pairs based on the samplemap DataFrame."""
        if 'condition' in df.columns:
            unique_conditions = df['condition'].dropna().unique()
            pairs = [(c1, c2) for c1, c2 in itertools.permutations(unique_conditions, 2)]
            self.cond_pairs = pairs
            pairs_str = [f"{c1}_VS_{c2}" for c1, c2 in pairs]
            self.condpairname_select.options = ["No conditions"] + pairs_str

    def _update_fc_visualizer(self):
        """Update FoldChangeVisualizer with current settings."""
        if hasattr(self, 'fc_visualizer') and self.cond1 and self.cond2:
            try:
                # Save DataFrame temporarily if needed
                if self.state.samplemap_df is not None and not self.state.samplemap_df.empty:
                    temp_dir = os.path.join(self.results_dir_input.value, 'temp')
                    os.makedirs(temp_dir, exist_ok=True)
                    temp_path = os.path.join(temp_dir, 'current_samplemap.tsv')
                    self.state.samplemap_df.to_csv(temp_path, sep='\t', index=False)

                    # Initialize visualizer with file path
                    self.fc_visualizer = aq_plot_fcviz.FoldChangeVisualizer(
                        condition1=self.cond1,
                        condition2=self.cond2,
                        results_directory=self.results_dir_input.value,
                        samplemap_file=temp_path,  # Use file path instead of DataFrame
                        tree_level=self.tree_level_select.value
                    )
                else:
                    # Initialize without samplemap if none available
                    self.fc_visualizer = aq_plot_fcviz.FoldChangeVisualizer(
                        condition1=self.cond1,
                        condition2=self.cond2,
                        results_directory=self.results_dir_input.value,
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
            self.condpairname_select.options = ["No conditions"] + pairs_str
        else:
            self.condpairname_select.options = ["No conditions"]

    def _on_condpair_selected(self, event):
        """Called when user selects a condition pair."""
        selected_str = event.new
        if selected_str == "No conditions":
            self.cond1 = None
            self.cond2 = None
            self._clear_plots()
            return

        if "_VS_" not in selected_str:
            return

        self.cond1, self.cond2 = selected_str.split("_VS_")
        self._update_data_for_condpair()
        self._build_volcano_plot()

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
                    height=600
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
            self.protein_input.value = event.new['points'][0]['text']

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
