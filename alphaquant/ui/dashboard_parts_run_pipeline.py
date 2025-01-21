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
import alphaquant.diffquant.diffutils as aqdiffutils
import alphaquant.run_pipeline as diffmgr
import alphaquant.ui.gui_textfields as gui_textfields

# If using Plotly in Panel
pn.extension('plotly')


class BaseWidget(param.Parameterized):
    """
    A base class to hold a common event trigger
    for Panel/Param watchers in child classes.
    """
    name = param.String()
    update_event = param.Integer(default=0, doc="Triggers re-computation")

    def __init__(self, name='', **params):
        super().__init__(name=name, **params)

    def trigger_dependency(self):
        """Increment update_event to force watchers to re-run."""
        self.update_event += 1


class HeaderWidget(param.Parameterized):
    """
    Header layout with project logos and links.
    """
    title = param.String()
    img_folder_path = param.String()
    github_url = param.String()

    def __init__(self, title, img_folder_path, github_url, **params):
        super().__init__(
            title=title,
            img_folder_path=img_folder_path,
            github_url=github_url,
            **params
        )
        self._create_panes()

    def _create_panes(self):
        """
        Initialize all Pane objects used in the header.
        """
        self.header_title = pn.pane.Markdown(
            f'# {self.title}',
            sizing_mode='stretch_width'
        )
        self.biochem_logo_path = os.path.join(self.img_folder_path, "mpi_logo.png")
        self.mpi_logo_path = os.path.join(self.img_folder_path, "max-planck-gesellschaft.jpg")
        self.github_logo_path = os.path.join(self.img_folder_path, "github.png")

        # Logos (png & jpg) with links
        self.mpi_biochem_logo = pn.pane.PNG(
            self.biochem_logo_path,
            link_url='https://www.biochem.mpg.de/mann',
            width=60, height=60
        )
        self.mpi_logo = pn.pane.JPG(
            self.mpi_logo_path,
            link_url='https://www.biochem.mpg.de/en',
            height=60, width=60
        )
        self.github_logo = pn.pane.PNG(
            self.github_logo_path,
            link_url=self.github_url,
            height=60
        )

    def create(self):
        """
        Return a layout (Panel Row) of the header widgets.
        """
        return pn.Row(
            self.mpi_biochem_logo,
            self.mpi_logo,
            self.header_title,
            self.github_logo,
            sizing_mode='stretch_width',
            margin=(5, 10, 5, 10)
        )


class MainWidget(param.Parameterized):
    """
    Create a layout for tool description and an optional manual download.
    """
    description = param.String()
    manual_path = param.String()

    def __init__(self, description, manual_path, **params):
        super().__init__(description=description, manual_path=manual_path, **params)
        self._create_widgets()

    def _create_widgets(self):
        """
        Initialize your markdown description and file downloader.
        """
        self.project_description = pn.pane.Markdown(
            self.description,
            margin=(10, 0, 10, 0),
            sizing_mode='stretch_width'
        )
        self.manual = pn.widgets.FileDownload(
            file=self.manual_path,
            label='Download Manual',
            button_type='default',
            auto=True,
            margin=(5, 5, 5, 5)
        )

    def create(self):
        """
        Return a simple column layout with the main description (and manual download).
        """
        layout = pn.Column(
            self.project_description,
            self.manual,
            sizing_mode='stretch_width',
            margin=(5, 5, 5, 5)
        )
        return layout


class RunPipeline(BaseWidget):
    """
    Widget to gather file inputs, define condition pairs, and run an analysis pipeline.
    Includes advanced configuration options.
    """
    def __init__(self, **params):
        super().__init__(**params)
        self._setup_matplotlib()
        self._make_widgets()
        self.layout = None

    def _make_widgets(self):
        # Original widgets
        self.path_analysis_file = pn.widgets.TextInput(
            name='Analysis file:',
            placeholder='Path to MQ/Spectronaut/DIA-NN file',
            width=700,
            sizing_mode='fixed'
        )
        self.path_output_folder = pn.widgets.TextInput(
            name='Output folder:',
            placeholder='Path to output folder',
            width=700,
            sizing_mode='fixed'
        )
        self.samplemap_title = pn.pane.Markdown('**Load an experiments-to-conditions file**')
        self.samplemap = pn.widgets.FileInput(
            accept='.tsv,.csv,.txt',
            margin=(5, 5, 10, 0)
        )
        self.samplemap_table = pn.widgets.Tabulator(
            layout='fit_data_fill',
            height=250,
            show_index=False,
            width=500,
            margin=(5, 5, 5, 5)
        )
        self.assign_cond_pairs = pn.widgets.CrossSelector(
            width=600,
            height=250,
            margin=(5, 5, 5, 5),
            name='Select condition pairs'
        )

        # New configuration widgets
        self.modification_type = pn.widgets.TextInput(
            name='Modification type:',
            placeholder='e.g., [Phospho (STY)] for Spectronaut',
            width=300
        )
        self.input_type = pn.widgets.TextInput(
            name='Input type:',
            placeholder='Type of quantitative information',
            width=300
        )
        self.organism = pn.widgets.TextInput(
            name='Organism:',
            placeholder='e.g., human, mouse',
            width=300
        )
        self.minrep_both = pn.widgets.IntInput(
            name='Min replicates (both conditions):',
            value=2,
            start=1,
            width=300
        )
        self.min_num_ions = pn.widgets.IntInput(
            name='Min number of ions per peptide:',
            value=1,
            start=1,
            width=300
        )
        self.minpep = pn.widgets.IntInput(
            name='Min peptides per protein:',
            value=1,
            start=1,
            width=300
        )
        self.cluster_threshold_pval = pn.widgets.FloatInput(
            name='Clustering p-value threshold:',
            value=0.001,
            start=0,
            end=1,
            width=300
        )
        self.volcano_fdr = pn.widgets.FloatInput(
            name='Volcano plot FDR:',
            value=0.05,
            start=0,
            end=1,
            width=300
        )
        self.volcano_fcthresh = pn.widgets.FloatInput(
            name='Volcano plot fold change threshold:',
            value=0.5,
            start=0,
            width=300
        )

        # Boolean switches
        self.switches = {
            'multicond_median_analysis': pn.widgets.Switch(name='Use median condition analysis', value=False),
            'use_ml': pn.widgets.Switch(name='Enable machine learning', value=True),
            'take_median_ion': pn.widgets.Switch(name='Use median-centered ions', value=True),
            'perform_ptm_mapping': pn.widgets.Switch(name='Enable PTM mapping', value=False),
            'perform_phospho_inference': pn.widgets.Switch(name='Enable phospho inference', value=False),
            'outlier_correction': pn.widgets.Switch(name='Enable outlier correction', value=True),
            'normalize': pn.widgets.Switch(name='Enable normalization', value=True),
            'use_iontree_if_possible': pn.widgets.Switch(name='Use ion tree when possible', value=True),
            'write_out_results_tree': pn.widgets.Switch(name='Write results tree', value=True),
            'use_multiprocessing': pn.widgets.Switch(name='Enable multiprocessing', value=False),
            'runtime_plots': pn.widgets.Switch(name='Generate runtime plots', value=True),
            'protnorm_peptides': pn.widgets.Switch(name='Enable protein-level normalization', value=True),
        }

        # Run pipeline widgets
        self.run_pipeline_button = pn.widgets.Button(
            name='Run pipeline',
            button_type='primary',
            height=35,
            width=170,
            margin=(10, 0, 0, 5)
        )
        self.run_pipeline_progress = pn.indicators.Progress(
            active=False,
            bar_color='light',
            width=170,
            margin=(5, 0, 5, 5)
        )
        self.visualize_data_button = pn.widgets.Button(
            name='Visualize data',
            button_type='success',
            height=35,
            width=170,
            margin=(10, 0, 0, 5)
        )
        self.run_pipeline_error = pn.pane.Alert(
            alert_type="danger",
            visible=False,
            margin=(5, 5, 5, 5),
        )

        # Attach watchers
        self.path_analysis_file.param.watch(
            self._activate_after_analysis_file_upload, 'value'
        )
        self.samplemap.param.watch(self._update_samplemap, 'value')
        self.samplemap_table.param.watch(
            self._add_conditions_for_assignment, 'value'
        )
        self.run_pipeline_button.param.watch(self._run_pipeline, 'clicks')
        self.visualize_data_button.param.watch(self._visualize_data, 'clicks')

    def create(self):
        """
        Build and return the main layout for the pipeline widget.
        """
        # Configuration card with all the settings
        config_card = pn.Card(
            pn.Column(
                "### Basic Settings",
                self.modification_type,
                self.input_type,
                self.organism,
                pn.layout.Divider(),
                "### Threshold Settings",
                self.minrep_both,
                self.min_num_ions,
                self.minpep,
                self.cluster_threshold_pval,
                self.volcano_fdr,
                self.volcano_fcthresh,
                pn.layout.Divider(),
                "### Analysis Options",
                pn.Column(*list(self.switches.values())),
            ),
            title='Advanced Configuration',
            collapsed=True,
            margin=(5, 5, 5, 5)
        )

        left_col = pn.Column(
            "### Input Files",
            self.path_analysis_file,
            self.path_output_folder,
            pn.Spacer(height=15),
            self.samplemap_title,
            self.samplemap,
            pn.Spacer(height=10),
            pn.Card(
                self.samplemap_table,
                title='Experiment → Condition map',
                collapsed=False,
                margin=(5, 5, 5, 5)
            ),
            pn.Card(
                self.assign_cond_pairs,
                title='Condition pairs for analysis',
                collapsed=True,
                margin=(5, 5, 5, 5)
            ),
            config_card,
            sizing_mode='stretch_width'
        )

        right_col = pn.Column(
            "### Instructions",
            gui_textfields.Descriptions.project_instruction,
            gui_textfields.Cards.spectronaut,
            gui_textfields.Cards.diann,
            gui_textfields.Cards.alphapept,
            gui_textfields.Cards.maxquant,
            pn.Spacer(height=10),
            pn.Row(self.run_pipeline_button, self.run_pipeline_progress),
            self.visualize_data_button,
            self.run_pipeline_error,
            sizing_mode='stretch_width'
        )

        self.layout = pn.Card(
            pn.Row(
                left_col,
                right_col,
                sizing_mode='stretch_width'
            ),
            title='Run Pipeline | Visualize Data',
            header_color='#333',
            header_background='#eaeaea',
            sizing_mode='stretch_width',
            margin=(10, 10, 10, 10)
        )
        return self.layout

    def _setup_matplotlib(self):
        """Configure matplotlib to use a non-GUI backend"""
        import matplotlib
        matplotlib.use('agg')  # Use the 'agg' backend which doesn't require GUI
        import matplotlib.pyplot as plt
        plt.ioff()  # Turn off interactive mode

    def _run_pipeline(self, *events):
        """Run the alphaquant pipeline when the button is clicked."""
        if not hasattr(self, 'data') or self.data.empty:
            self.run_pipeline_error.object = "No valid data loaded."
            self.run_pipeline_error.visible = True
            return

        self.run_pipeline_progress.active = True
        self.run_pipeline_error.visible = False

        # Get condition combinations from the CrossSelector
        if self.assign_cond_pairs.value:
            cond_combinations = [
                tuple(pair.split('_vs_'))
                for pair in self.assign_cond_pairs.value
            ]
        else:
            cond_combinations = [
                tuple(pair.split('_vs_'))
                for pair in self.assign_cond_pairs.options
            ]

        try:
            # Collect all configuration parameters
            pipeline_params = {
                'input_file': self.path_analysis_file.value,
                'samplemap_df': self.samplemap_table.value,
                'results_dir': self.path_output_folder.value,
                'condpairs_list': cond_combinations,
                'modification_type': self.modification_type.value or None,
                'input_type_to_use': self.input_type.value or None,
                'organism': self.organism.value or None,
                'minrep_both': self.minrep_both.value,
                'min_num_ions': self.min_num_ions.value,
                'minpep': self.minpep.value,
                'cluster_threshold_pval': self.cluster_threshold_pval.value,
                'volcano_fdr': self.volcano_fdr.value,
                'volcano_fcthresh': self.volcano_fcthresh.value,
            }

            # Add all switch values
            pipeline_params.update({
                key: switch.value for key, switch in self.switches.items()
            })

            # Run the pipeline with all parameters
            diffmgr.run_pipeline(**pipeline_params)

        except Exception as e:
            self.run_pipeline_error.object = f"Error running pipeline: {e}"
            self.run_pipeline_error.visible = True

        self.trigger_dependency()
        self.run_pipeline_progress.active = False

    # Other methods remain the same as in your original code
    def _activate_after_analysis_file_upload(self, *events):
        self._set_default_output_folder()
        self._import_exp_data()
        self._extract_sample_names()

    def _set_default_output_folder(self):
        if (not self.path_output_folder.value) and self.path_analysis_file.value:
            base_path = os.path.dirname(self.path_analysis_file.value)
            self.path_output_folder.value = os.path.join(base_path, 'results')

    def _import_exp_data(self):
        if self.path_analysis_file.value:
            try:
                self.data = aqdiffutils.import_data(input_file=self.path_analysis_file.value)
            except Exception as e:
                self.run_pipeline_error.object = f"Error importing data: {e}"
                self.run_pipeline_error.visible = True
                self.data = pd.DataFrame()
        else:
            self.data = pd.DataFrame()

    def _extract_sample_names(self):
        if hasattr(self, 'data') and not self.data.empty:
            sample_names = aqdiffutils.get_samplenames_from_input_df(self.data)
            sorted_names = self.natural_sort(sample_names)
            self.samplemap_table.value = pd.DataFrame({
                'sample': sorted_names,
                'condition': [''] * len(sorted_names)
            })

    def _update_samplemap(self, *events):
        if not self.samplemap.value:
            return
        file_ext = os.path.splitext(self.samplemap.filename)[-1].lower()
        sep = ',' if file_ext == '.csv' else '\t'
        try:
            df = pd.read_csv(
                StringIO(self.samplemap.value.decode('utf-8')),
                sep=sep,
                dtype=str
            )
            self.samplemap_table.value = df
        except Exception as e:
            self.run_pipeline_error.object = f"Error reading sample map: {e}"
            self.run_pipeline_error.visible = True

    def _add_conditions_for_assignment(self, *events):
        if self.samplemap_table.value is None:
            return
        df = self.samplemap_table.value
        if 'condition' in df.columns:
            unique_condit = df['condition'].dropna().unique()
            comb_condit = [
                '_vs_'.join(comb)
                for comb in itertools.permutations(unique_condit, 2)
            ]
            self.assign_cond_pairs.options = comb_condit

    def _visualize_data(self, *events):
        self.run_pipeline_error.visible = False
        self.trigger_dependency()

    def natural_sort(self, l):
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(l, key=alphanum_key)


class Tabs(param.Parameterized):
    """
    Manages switching from the pipeline stage to one or more analysis tabs.
    """
    pipeline = param.ClassSelector(class_=RunPipeline)
    main_tabs = param.ClassSelector(class_=pn.Tabs, allow_None=True)
    current_layout = param.Parameter()

    def __init__(self, pipeline, **params):
        super().__init__(pipeline=pipeline, **params)
        self.main_tabs = None
        self.current_layout = None
        # Watch pipeline.update_event to create the tabs after pipeline has run
        self.pipeline.param.watch(self._create_tabs_if_possible, 'update_event')

    def _create_tabs_if_possible(self, *events):
        """
        If the user has clicked 'visualize data', show the analysis tabs.
        """
        if (self.pipeline.path_output_folder.value
                and self.pipeline.visualize_data_button.clicks > 0):
            try:
                self._build_tabs()
                # Update the current layout
                self.current_layout = self.main_tabs
            except Exception as e:
                print(f"Error creating visualization tabs: {str(e)}")
                self.current_layout = pn.pane.Markdown(
                    f"Error creating visualization: {str(e)}"
                )

    def _build_tabs(self):
        """
        Build a Panel Tabs layout with visualization components.
        """
        if self.main_tabs is None:
            self.main_tabs = pn.Tabs(
                tabs_location='above',
                sizing_mode='stretch_width',
                margin=(10, 10, 10, 10)
            )

            # Initialize visualization components
            try:
                # Single Comparison tab
                single_comp = SingleComparison(
                    self.pipeline.path_output_folder.value,
                    self.pipeline.samplemap_table.value
                )
                self.main_tabs.append(
                    ('Single Comparison', single_comp.layout)
                )
                # Could add more tabs here for multiple comparisons, etc.
            except Exception as e:
                print(f"Error initializing comparison tab: {str(e)}")
                self.main_tabs.append(
                    ('Error', pn.pane.Markdown(f"Error creating visualization: {str(e)}"))
                )

    def create(self):
        """
        Return the current view - either tabs or a message.
        """
        if self.current_layout is not None:
            return self.current_layout
        else:
            return pn.pane.Markdown(
                "## No results yet.\nPlease run the pipeline and click **Visualize data**.",
                sizing_mode='stretch_width'
            )

# ----------------
# BUILD DASHBOARD
# ----------------
def build_dashboard():
    """
    Example function to build the overall dashboard layout in a FastListTemplate.
    """
    header = HeaderWidget(
        title="AlphaQuant Dashboard",
        img_folder_path="./assets",
        github_url="https://github.com/<my_repo>"
    )
    main_text = MainWidget(
        description=(
            "Welcome to our analysis dashboard. "
            "Please load your data and run the pipeline."
        ),
        manual_path="path/to/manual.pdf"
    )
    pipeline = RunPipeline()
    tab_manager = Tabs(pipeline)

    template = pn.template.FastListTemplate(
        title="AlphaQuant Analysis",
        sidebar=[   # If you want a sidebar, you can put items here
            # "## Sidebar Title",
            # pn.widgets.Select(options=['Item1','Item2']),
        ],
        main=[
            header.create(),
            pn.layout.Divider(),
            main_text.create(),
            pipeline.create(),
            tab_manager.create()
        ],
        theme='dark',          # or 'default'
        main_max_width="1200px"
    )
    return template

# If run as a script, you can do:
# if __name__ == "__main__":
#     dash = build_dashboard()
#     dash.servable()
