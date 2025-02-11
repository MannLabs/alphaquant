import param
import panel as pn
import pandas as pd
import matplotlib.pyplot as plt
import os
from io import StringIO
import itertools
import tempfile
import glob
import re

import alphaquant.plotting.fcviz as aq_plot_fcviz
import alphaquant.plotting.alphamapviz as aq_plot_proteoform
import alphaquant.utils.proteoform_utils as aq_proteoform_utils

class ProteoformPlottingTab(param.Parameterized):
    """
    A param-based Panel class for proteoform visualization.
    """
    # Public Parameters
    results_dir = param.String(doc="Path to the results folder containing *_VS_*.results.tsv files")
    samplemap_file = param.String(doc="Path to the samplemap file")

    # Dynamic widgets
    condpairname_select = param.ClassSelector(class_=pn.widgets.Select)
    protein_input = param.ClassSelector(class_=pn.widgets.AutocompleteInput)
    proteoform_view_select = param.ClassSelector(class_=pn.widgets.Select)

    # Add new parameter for the selected proteoform
    selected_proteoform = param.String(default='')

    # Add new widgets for organism and protein identifier
    organism_select = pn.widgets.Select(
        name="Organism",
        options=['Human', 'Mouse', 'Yeast'],
        value='Human',
        width=200
    )

    protein_id_select = pn.widgets.Select(
        name="Protein Identifier",
        options=['gene_symbol', 'uniprot_id'],
        value='gene_symbol',
        width=200
    )

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
            placeholder="Type protein name or click a row in the table above to visualize peptide fold changes...",
            min_characters=2,
            disabled=True,
            width=400,
            margin=(5, 5, 5, 5)
        )
        self.protein_input.param.watch(self._on_protein_selected, 'value')

        self.proteoform_view_select = pn.widgets.Select(
            name="Visualization Type",
            options=['Sequence Plot'],
            value='Sequence Plot',
            width=200
        )

        # Input fields for paths
        self.results_dir_input = pn.widgets.TextInput(
            name='Results Directory:',
            value=self.results_dir,
            placeholder='Enter path to results directory',
            width=600
        )
        self.results_dir_input.param.watch(self.on_results_dir_changed, 'value')

        # Replace file upload with text input for samplemap
        self.samplemap_input = pn.widgets.TextInput(
            name='Sample Map File:',
            value=self.samplemap_file,
            placeholder='Enter path to sample map file',
            width=600
        )
        self.samplemap_input.param.watch(self._on_samplemap_changed, 'value')

        # Add new widget for proteoform table
        self.proteoform_table = pn.widgets.Tabulator(
            pagination='remote',
            page_size=10,
            sizing_mode='stretch_width',
            height=300,
            selectable=1,
            selection=[],
            show_index=False,
        )
        self.proteoform_table.on_click(self._on_proteoform_selected)

        # Plot panes
        self.proteoform_plot_pane = pn.Column(
            pn.pane.Markdown("### Visualization will appear here when you select a protein"),
            sizing_mode='stretch_width'
        )

        # Construct layout
        self.main_layout = pn.Column(
            "## Outlier Peptide Visualization",
            self.results_dir_input,
            self.samplemap_input,
            pn.Row(self.organism_select, self.protein_id_select),  # Add new widgets
            self.condpairname_select,
            self.proteoform_table,
            pn.Row(self.proteoform_view_select),
            self.protein_input,
            self.proteoform_plot_pane,
            sizing_mode='stretch_width'
        )

        # Initialize if directory provided
        self._extract_condpairs()

    def panel(self):
        """Return the main panel layout."""
        return self.main_layout

    def on_results_dir_changed(self, event):
        """Handle changes to results directory."""
        if event.new:
            self.results_dir = event.new
            self._extract_condpairs()

    def _extract_condpairs(self):
        """Look for '*_VS_*.proteoforms.tsv' in the results_dir and update the condition pairs."""
        if not self.results_dir or not os.path.isdir(self.results_dir):
            self.condpairname_select.options = ["No conditions"]
            return

        pattern = os.path.join(self.results_dir, "*_VS_*.proteoforms.tsv")
        files = glob.glob(pattern)

        cond_pairs = []
        for f in files:
            basename = os.path.basename(f)
            match = re.match(r'(.*?)_VS_(.*?)\.proteoforms\.tsv$', basename)
            if match:
                cond1, cond2 = match.group(1), match.group(2)
                cond_pairs.append((cond1, cond2))

        if cond_pairs:
            pairs_str = [f"{c1}_VS_{c2}" for c1, c2 in cond_pairs]
            self.condpairname_select.options = ["No conditions"] + pairs_str
        else:
            self.condpairname_select.options = ["No conditions"]

    def _on_condpair_selected(self, event):
        """Handle condition pair selection."""
        if event.new and event.new != "No conditions":
            condition1, condition2 = event.new.split('_VS_')
            results_file = os.path.join(
                self.results_dir,
                f"{condition1}_VS_{condition2}.proteoforms.tsv"
            )

            try:
                # Load and filter proteoforms
                proteoforms_df = pd.read_csv(results_file, sep='\t')
                filtered_df = aq_proteoform_utils.filter_proteoform_df(proteoforms_df)

                # Drop specified columns
                columns_to_drop = ['is_reference', 'peptides', 'log2fc',
                                 'proteoform_pval', 'proteoform_fcfc', 'fcdiff',
                                 'proteoform_fdr']
                filtered_df = filtered_df.drop(columns=[col for col in columns_to_drop if col in filtered_df.columns])

                # Update table
                self.proteoform_table.value = filtered_df
                self.proteoform_table.visible = True

                # Update protein input options
                protein_ids = filtered_df['protein'].unique().tolist()
                self.protein_input.options = protein_ids
                self.protein_input.disabled = False

                print("Selected condition pair:", event.new)
                print(f"Parsed conditions: {condition1=}, {condition2=}")

                print("Updating visualizers")
                # Initialize both visualizers
                self.amap_visualizer = aq_plot_proteoform.AlphaMapVisualizer(
                    condition1=condition1,
                    condition2=condition2,
                    results_directory=self.results_dir,
                    samplemap_file=self.samplemap_file,
                    protein_identifier=self.protein_id_select.value,
                    organism=self.organism_select.value
                )

                self.fc_visualizer = aq_plot_fcviz.FoldChangeVisualizer(
                    condition1=condition1,
                    condition2=condition2,
                    results_directory=self.results_dir,
                    samplemap_file=self.samplemap_file,
                    organism=self.organism_select.value,
                    protein_identifier=self.protein_id_select.value,
                    order_along_protein_sequence=True,
                    figsize=(6, 4)  # Smaller figure size for fold change plot
                )

            except Exception as e:
                print("Error occurred:", str(e))
                print("Exception type:", type(e))
                import traceback
                print("Traceback:", traceback.format_exc())
                self.protein_input.disabled = True
                self.protein_input.options = []
                error_msg = f"Error loading proteoforms file: {str(e)}"
                self.proteoform_plot_pane.clear()
                self.proteoform_plot_pane.append(pn.pane.Markdown(f"### Error\n{error_msg}"))

    def _on_protein_selected(self, event):
        """Handle protein selection."""
        if event.new:
            # Clear existing plots
            self.proteoform_plot_pane.clear()

            # Get sequence plot from AlphaMapVisualizer
            _, alphamap_go_fig = self.amap_visualizer.visualize_protein(event.new)

            # Get fold change plot from FoldChangeVisualizer
            fc_fig = self.fc_visualizer.plot_protein(event.new)

            # Add both plots to the pane
            if fc_fig:
                self.proteoform_plot_pane.append(pn.pane.Matplotlib(fc_fig, tight=True))
            if alphamap_go_fig:
                self.proteoform_plot_pane.append(pn.pane.Plotly(alphamap_go_fig))

    def _on_samplemap_changed(self, event):
        """Handle changes to samplemap file path."""
        if event.new:
            self.samplemap_file = event.new
            try:
                # Verify the file exists and can be read
                df = pd.read_csv(self.samplemap_file, sep='\t', dtype=str)
                num_samples = len(df)
                num_conditions = len(df['condition'].unique()) if 'condition' in df.columns else 0
                print(f"Loaded sample map with {num_samples} samples and {num_conditions} conditions")
            except Exception as e:
                print(f"Error loading sample map: {str(e)}")

    def on_samplemap_df_changed(self, new_df):
        """Handle changes to samplemap DataFrame from other components."""
        if not new_df.empty:
            # Update status
            num_samples = len(new_df)
            num_conditions = len(new_df['condition'].unique()) if 'condition' in new_df.columns else 0
            status_text = f"Sample Map: Loaded {num_samples} samples, {num_conditions} conditions"
            self.samplemap_input.value = status_text

            # Update condition pairs and other visualizations
            self._update_condition_pairs_from_df(new_df)
        else:
            self.samplemap_input.value = "Sample Map: No sample map loaded"

    def _update_condition_pairs_from_df(self, df):
        """Update condition pairs based on the samplemap DataFrame."""
        if 'condition' in df.columns:
            unique_conditions = df['condition'].dropna().unique()
            pairs = [(c1, c2) for c1, c2 in itertools.permutations(unique_conditions, 2)]
            pairs_str = [f"{c1}_VS_{c2}" for c1, c2 in pairs]
            self.condpairname_select.options = ["No conditions"] + pairs_str

    def _load_protein_identifiers(self, results_file):
        """Load protein identifiers from results file and update the protein input widget."""
        try:
            proteoforms_df = pd.read_csv(results_file, sep='\t')
            protein_ids = sorted(proteoforms_df['gene_symbol'].unique().tolist())
            self.protein_input.options = protein_ids
            self.protein_input.disabled = False
            return protein_ids
        except Exception as e:
            self.protein_input.disabled = True
            self.protein_input.options = []
            raise Exception(f"Failed to load protein identifiers: {str(e)}")

    def _on_proteoform_selected(self, event):
        """Handle proteoform selection from table."""
        if hasattr(event, 'row'):
            row_data = self.proteoform_table.value.iloc[event.row]
            selected_protein = row_data.get('protein')
            if selected_protein:
                self.protein_input.value = selected_protein
                # Directly update the plot without requiring click on protein input
                self._on_protein_selected(param.Event(type='selection', new=selected_protein))
