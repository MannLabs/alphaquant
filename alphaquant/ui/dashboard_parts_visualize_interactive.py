import os
import re
from io import StringIO
import itertools

import param
import panel as pn
import pandas as pd
import numpy as np
from scipy import stats
import anytree
import matplotlib
matplotlib.use('agg')

# alphaquant imports
import alphaquant.utils.utils as aqutils
import alphaquant.plotting.base_functions as aqplot
import alphaquant.config.variables as aq_variables

# If using Plotly in Panel
pn.extension('plotly')


class SingleComparison(param.Parameterized):
    """
    Creates a panel with single comparison-based plots (volcano, beeswarm, etc.).
    """
    output_folder = param.String()
    sample_to_cond = param.DataFrame()

    condpairs_selector = param.ClassSelector(class_=pn.widgets.Select)
    protein = param.ClassSelector(class_=pn.widgets.AutocompleteInput)

    def __init__(self, output_folder, sample_to_cond, **params):
        super().__init__(output_folder=output_folder, sample_to_cond=sample_to_cond, **params)
        self._initialize_widgets()
        self._init_data()
        self.layout = self.create()

    def _initialize_widgets(self):
        self.condpairs_selector = pn.widgets.Select(
            name='Select a pair of conditions:',
            options=['No conditions'],
            width=400,
            margin=(5, 5, 5, 5)
        )
        self.protein = pn.widgets.AutocompleteInput(
            name='Select protein:',
            placeholder="Type the first letters",
            min_characters=2,
            disabled=True,
            width=400,
            margin=(5, 5, 5, 5)
        )

        # Watchers
        self.condpairs_selector.param.watch(self._run_after_pair_cond_selection, 'value')
        self.protein.param.watch(self._visualize_after_protein_selection, 'value')

        self.volcano_plot = None
        self.result_df = pd.DataFrame()
        self.normalized_intensity_df = pd.DataFrame()
        self.iontree_condpair = None

    def _init_data(self):
        """
        Scan the output folder for .results.tsv files and populate the condpairs_selector.
        """
        if not os.path.isdir(self.output_folder):
            return
        results_files = [
            f.replace(".results.tsv", "").replace('VS', 'vs')
            for f in os.listdir(self.output_folder)
            if f.endswith(".results.tsv")
        ]
        self.condpairs_selector.options = ['No conditions'] + results_files

    def create(self):
        """
        Construct the SingleComparison layout.
        """
        layout = pn.Column(
            "## Single Comparison Plots",
            pn.Row(self.condpairs_selector, self.protein),
            pn.Spacer(height=10),
            # distribution plots (card)
            None,
            # volcano plot
            None,
            pn.layout.Divider(),
            # beeswarm & fold-change
            pn.Row(None, None),
            sizing_mode='stretch_width',
            margin=(10, 10, 10, 10)
        )
        return layout

    def _run_after_pair_cond_selection(self, *events):
        layout = self.layout
        if self.condpairs_selector.value == 'No conditions':
            layout[2] = None  # distribution plots
            layout[3] = None  # volcano
            layout[5][0] = None  # beeswarm
            layout[5][1] = None  # fold-change
            return

        # Parse the condition pair
        self.cond1, self.cond2 = self.condpairs_selector.value.split(aq_variables.CONDITION_PAIR_SEPARATOR)
        self.iontree_condpair = aqutils.read_condpair_tree(self.cond1, self.cond2, self.output_folder)

        # Read result & normalized data
        self.result_df = aqplot.get_diffresult_dataframe(
            self.cond1, self.cond2, results_folder=self.output_folder
        )
        self.normalized_intensity_df = aqplot.get_normed_peptides_dataframe(
            self.cond1, self.cond2, results_folder=self.output_folder
        )

        # Initialize the helper to get protein-based intensities
        quantinfo = aqplot.CondpairQuantificationInfo().init_from_loaded_tables(
            diffresults_df=self.result_df,
            normed_df=self.normalized_intensity_df,
            condpair_root_node=self.iontree_condpair,
            samplemap_df=self.sample_to_cond
        )
        self.protein_df_getter = aqplot.ProteinIntensityDataFrameGetter(quantinfo)

        # Populate protein autocomplete
        self.protein.options = self.result_df.protein.dropna().unique().tolist()
        self.protein.disabled = False

        # Some distribution plots
        c1_normed = aqplot.subset_normed_peptides_df_to_condition(
            self.cond1, self.sample_to_cond, self.normalized_intensity_df
        )
        c2_normed = aqplot.subset_normed_peptides_df_to_condition(
            self.cond2, self.sample_to_cond, self.normalized_intensity_df
        )

        dist_plots_card = pn.Card(
            "# Log2 Fold Change Distributions",
            pn.pane.Plotly(
                aqplot.plot_withincond_fcs_plotly(
                    self.normalized_intensity_df,
                    'All samples'
                ),
                config={'responsive': True}
            ),
            pn.pane.Plotly(
                aqplot.plot_withincond_fcs_plotly(
                    c1_normed,
                    f'Condition: {self.cond1}'
                ),
                config={'responsive': True}
            ),
            pn.pane.Plotly(
                aqplot.plot_withincond_fcs_plotly(
                    c2_normed,
                    f'Condition: {self.cond2}'
                ),
                config={'responsive': True}
            ),
            pn.pane.Plotly(
                aqplot.plot_betweencond_fcs_plotly(
                    c1_normed, c2_normed,
                    'Between Conditions'
                ),
                config={'responsive': True}
            ),
            collapsed=True,
            sizing_mode='stretch_width',
            margin=(5, 5, 5, 5)
        )

        self.volcano_plot = pn.pane.Plotly(
            aqplot.plot_volcano_plotly(self.result_df),
            config={'responsive': True}
        )
        self.volcano_plot.param.watch(self._on_volcano_click, 'click_data')

        layout[2] = dist_plots_card
        layout[3] = self.volcano_plot
        layout[5][0] = None
        layout[5][1] = None

    def _on_volcano_click(self, event):
        """
        When user clicks a point in the volcano, update the protein selector.
        """
        if event and 'points' in event.obj and event.obj['points']:
            self.protein.value = event.obj['points'][0]['text']

    def _visualize_after_protein_selection(self, *events):
        layout = self.layout
        if not self.protein.value or self.protein.value not in self.result_df.protein.values:
            layout[5][0] = None
            layout[5][1] = None
            return

        # Melted DF for chosen protein
        melted_df = self.protein_df_getter.get_melted_protein_ion_intensity_table(self.protein.value)
        self.protein_df = self.protein_df_getter.get_protein_diffresults(self.protein.value)

        # Between-cond table
        fc_df = aqplot.get_betweencond_fcs_table(melted_df, self.cond1, self.cond2)

        # Possibly find the protein node in the condpair tree
        protnode = None
        if self.iontree_condpair is not None:
            protnode = anytree.find(
                self.iontree_condpair,
                filter_=lambda x: x.name == self.protein.value,
                maxlevel=2
            )

        beeswarm_pane = pn.pane.Plotly(
            aqplot.beeswarm_ion_plot_plotly(melted_df, self.protein_df),
            config={'responsive': True}
        )
        fc_pane = pn.pane.Plotly(
            aqplot.foldchange_ion_plot_plotly(fc_df, self.protein_df, protein_node=protnode),
            config={'responsive': True}
        )

        # Place them in the layout
        layout[5][0] = beeswarm_pane
        layout[5][1] = fc_pane


