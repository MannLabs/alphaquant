import os
import re
from io import StringIO
from itertools import permutations
import pandas as pd

# alphaquant important
import alphaquant.diffquant_utils as aqutils
import alphaquant.diff_analysis_manager as diffmgr
import alphaquant.visualizations as aqplot

# visualization
import panel as pn
import plotly.graph_objs as go

UPDATED = pn.widgets.IntInput(value=0)

class HeaderWidget(object):
    """This class creates a layout for the header of the dashboard with the name of the tool and all links to the MPI website, the MPI Mann Lab page and the GitHub repo.

    Parameters
    ----------
    title : str
        The name of the tool.

    Attributes
    ----------
    header_title : pn.pane.Markdown
        A Panel Markdown pane that returns the title of the tool.
    mpi_biochem_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the MPI Biochemisty logo and makes the image clickable with the link to the official website.
    mpi_logo : pn.pane.JPG
        A Panel JPG pane that embeds a jpg image file of the MPI Biochemisty logo and makes the image clickable with the link to the official website.
    github_logo : pn.pane.PNG
        A Panel PNG pane that embeds a png image file of the GitHub logo and makes the image clickable with the link to the GitHub repository of the project.

    """

    def __init__(self,
        title,
        img_folder_path,
        github_url
        ):
        self.header_title = pn.pane.Markdown(
            f'# {title}',
            sizing_mode='stretch_width',
        )
        self.biochem_logo_path = os.path.join(
            img_folder_path,
            "mpi_logo.png"
        )
        self.mpi_logo_path = os.path.join(
            img_folder_path,
            "max-planck-gesellschaft.jpg"
        )
        self.github_logo_path = os.path.join(
            img_folder_path,
            "github.png"
        )

        self.mpi_biochem_logo = pn.pane.PNG(
            self.biochem_logo_path,
            link_url='https://www.biochem.mpg.de/en',
            width=60,
            height=60,
            align='start'
        )
        self.mpi_logo = pn.pane.JPG(
            self.mpi_logo_path,
            link_url='https://www.biochem.mpg.de/en',
            height=62,
            embed=True,
            width=62,
            margin=(5, 0, 0, 5),
            css_classes=['opt']
        )
        self.github_logo = pn.pane.PNG(
            self.github_logo_path,
            link_url=github_url,
            height=70,
            align='end'
        )


    def create(self):
        return pn.Row(
            self.mpi_biochem_logo,
            self.mpi_logo,
            self.header_title,
            self.github_logo,
            height=73,
            sizing_mode='stretch_width'
        )


class MainWidget(object):

    def __init__(self,
        description,
        manual_path
    ):
        self.project_description = pn.pane.Markdown(
            description,
            margin=(10, 0, 10, 0),
            css_classes=['main-part'],
            align='start'
        )
        self.manual = pn.widgets.FileDownload(
            file=manual_path,
            label='Download Manual',
            button_type='default',
            align='center',
            auto=True,
            height=31,
            width=200,
            margin=(0, 20, 0, 0)
        )


    def create(self):
        LAYOUT = pn.Row(
            self.project_description,
            pn.layout.HSpacer(width=500),
            self.manual,
            background='#eaeaea',
            align='center',
            sizing_mode='stretch_width',
            height=190,
            margin=(10, 8, 10, 8),
            css_classes=['background']
        )
        return LAYOUT


class RunAnalysis(object):

    def __init__(self):
        # DATA FILES
        self.path_analysis_file = pn.widgets.TextInput(
            name='Specify an analysis file:',
            placeholder='Enter the whole path to the MQ | Spectronaut | DIA-NN output file',
            width=900,
            sizing_mode='stretch_width',
            margin=(5, 15, 0, 15)
        )
        self.path_output_folder = pn.widgets.TextInput(
            name='Specify a path to the output folder:',
            placeholder='Enter the whole path to the output folder',
            width=900,
            sizing_mode='stretch_width',
            margin=(15, 15, 0, 15)
        )
        self.samplemap_title = pn.pane.Markdown(
            'Load an experiments-to-conditions file:',
            margin=(10,0,0,15)
        )
        self.samplemap = pn.widgets.FileInput(
            accept='.tsv,.csv,.txt',
            margin=(-10,0,5,15)
        )
        self.samplemap_table = pn.widgets.Tabulator(
            layout='fit_data_fill',
            height=300,
            show_index=False,
            width=570,
            align='center',
            margin=(15, 12, 10, 18)
        )
        self.assign_cond_pairs = pn.widgets.CrossSelector(
            width=870,
            height=300,
            align='center',
            margin=(15, 15, 15, 15)
        )
        # RUN PIPELINE
        self.run_pipeline_button = pn.widgets.Button(
            name='Run pipeline',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(0, 0, 0, 0)
        )
        self.run_pipeline_progress = pn.indicators.Progress(
            active=False,
            bar_color='light',
            width=250,
            align='center',
            margin=(-10, 0, 20, 0)
        )
        self.visualize_data_button = pn.widgets.Button(
            name='Visualize data',
            button_type='primary',
            height=31,
            width=250,
            align='center',
            margin=(40, 0, 0, 0)
        )
        self.run_pipeline_error = pn.pane.Alert(
            width=600,
            alert_type="danger",
            # object='test warning message',
            margin=(-20, 10, -5, 16),
        )
        # self.updated = pn.widgets.IntInput(value=0)
        self.data = None



    def create(self):
        LAYOUT = pn.Card(
            pn.Row(
                pn.Column(
                    self.path_analysis_file,
                    self.run_pipeline_error,
                    self.path_output_folder,
                    pn.Row(
                        pn.Column(
                            self.samplemap_title,
                            self.samplemap
                        ),
                        pn.Card(
                            self.samplemap_table,
                            header='Assign experiments to conditions manually',
                            collapsed=True,
                            margin=(20, 0, 20, 0),
                            width=601,
                        )
                    ),
                    pn.Card(
                        self.assign_cond_pairs,
                        header='Select condition pairs for the analysis',
                        collapsed=True,
                        margin=(5, 0, 20, 15),
                        width=901,
                    ),
                    margin=(20, 30, 10, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                pn.Column(
                    self.run_pipeline_button,
                    self.run_pipeline_progress,
                    self.visualize_data_button,
                    align='center',
                    margin=(100, 40, 0, 0),
                )
            ),
            title='Run Pipeline | Visualize data',
            collapsed=False,
            header_background='#eaeaea',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            height=1000,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        self.activate_after_analysis_file_upload = pn.depends(
            self.path_analysis_file.param.value,
            watch=True
        )(self.activate_after_analysis_file_upload)

        self.update_samplemap = pn.depends(
            self.samplemap.param.value,
            watch=True
        )(self.update_samplemap)

        self.add_conditions_for_assignment = pn.depends(
            self.samplemap_table.param.value,
            watch=True
        )(self.add_conditions_for_assignment)

        self.run_pipeline = pn.depends(
            self.run_pipeline_button.param.clicks,
            watch=True
        )(self.run_pipeline)

        self.visualize_data = pn.depends(
            self.visualize_data_button.param.clicks,
            watch=True
        )(self.visualize_data)

        return LAYOUT


    def activate_after_analysis_file_upload(self, *args):
        self.set_default_output_folder()
        self.import_exp_data()
        self.extract_sample_names()


    def natural_sort(self, l):
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(l, key = alphanum_key)


    def set_default_output_folder(self):
        if not self.path_output_folder.value:
            self.path_output_folder.value = os.path.join(
                os.path.dirname(self.path_analysis_file.value),
                'results'
            )


    def import_exp_data(self):
        self.data = aqutils.import_data(
            input_file=self.path_analysis_file.value,
            results_folder=self.path_output_folder.value
        )


    def extract_sample_names(self):
        sample_names = aqutils.get_samplenames(self.data)
        self.samplemap_table.value = pd.DataFrame(
            data={'sample': self.natural_sort(sample_names), 'condition': str()}
        )


    def update_samplemap(self, *args):
        file_ext = os.path.splitext(self.samplemap.filename)[-1]

        if file_ext=='.csv':
            sep=','
        else:
            sep='\t'

        self.samplemap_table.value = pd.read_csv(
            StringIO(str(self.samplemap.value, "utf-8")),
            sep=sep
        )


    def add_conditions_for_assignment(self, *args):
        unique_condit = self.samplemap_table.value.condition.unique()
        comb_condit = ['_vs_'.join(comb) for comb in permutations(unique_condit, 2)]
        self.assign_cond_pairs.options = comb_condit


    def run_pipeline(self, *args):
        global UPDATED
        self.run_pipeline_progress.active = True

        data_processed, samplemap_df_processed = aqutils.prepare_loaded_tables(
            self.data,
            self.samplemap_table.value
        )

        if self.assign_cond_pairs.options:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.assign_cond_pairs.value]
        else:
            cond_combinations = None

        diffmgr.run_pipeline(
            unnormed_df=data_processed,
            labelmap_df=samplemap_df_processed,
            results_dir=self.path_output_folder.value,
            condpair_combinations=cond_combinations
        )

        self.run_pipeline_progress.active = False
        UPDATED.value += 1


    def visualize_data(self, *args):
        global UPDATED
        print(UPDATED)
        UPDATED.value += 1


class Tabs(object):

    def __init__(self):
        self.LAYOUT = None


    def create(self, *args):
        UPDATED.param.watch(
            self.fill_layout(*args),
            'value'
        )
        return self.LAYOUT


    def fill_layout(self, *args):
        self.LAYOUT = pn.Tabs(
            *args,
            active=0,
            tabs_location='above',
            sizing_mode='stretch_width',
            margin=(5, 8, 10, 8)
        )


class MultipleComparison(object):

    def __init__(self, output_folder):
        self.condpairs_to_compare = pn.widgets.CrossSelector(
            options=['c2_vs_c12', 'c12_vs_c2'],
            width=870,
            height=300,
            align='center',
            margin=(15, 15, 15, 15)
        )
        self.output_folder = output_folder
        self.LAYOUT = None


    def create(self):
        self.condpairs_to_compare.options = self.extract_conditions_from_folder()
        self.condpairs_to_compare.param.watch(
            self.return_clustered_heatmap,
            'value'
        )
        self.LAYOUT = pn.Column(
            self.condpairs_to_compare,
        )
        return self.LAYOUT


    def extract_conditions_from_folder(self):
        return [f.replace(".results.tsv", "").replace('VS', 'vs') for f in os.listdir(self.output_folder) if re.match(r'.*results.tsv', f)]


    def return_clustered_heatmap(self):

        if self.condpairs_to_compare.value:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.condpairs_to_compare.value]
        else:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.condpairs_to_compare.options]

        overview_dataframe = aqplot.get_sample_overview_dataframe(
            results_folder=self.output_folder,
            condpairs_to_compare=cond_combinations
        )

    def plot_heatmap(df, title, colormap):
        """
        This function plots a simple Heatmap.
        """
        fig = go.Figure()
        fig.update_layout(
            title={
                'text': title,
                'y': 0.95,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            height=1000,
            width=700,
            annotations=[dict(xref='paper', yref='paper', showarrow=False, text='')],
            template='plotly_white'
        )
        fig.add_trace(
            go.Heatmap(
                z=df.values.tolist(),
                x = list(df.columns),
                y = list(df.index),
                colorscale=colormap
            )
        )
        return fig
