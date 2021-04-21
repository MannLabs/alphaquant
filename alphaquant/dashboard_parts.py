import os
import re
from io import StringIO
from itertools import permutations
import pandas as pd

# visualization
import panel as pn


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
        self.predefined_exp_to_cond_title = pn.pane.Markdown(
            'Load an experiments-to-conditions file:',
            margin=(10,0,0,15)
        )
        self.predefined_exp_to_cond = pn.widgets.FileInput(
            accept='.tsv,.csv',
            margin=(-10,0,5,15)
        )
        self.df_exp_to_cond = pn.widgets.Tabulator(
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
        self.updated = pn.widgets.IntInput(value=0)


    def create(self):
        LAYOUT = pn.Card(
            pn.Row(
                pn.Column(
                    self.path_analysis_file,
                    self.run_pipeline_error,
                    self.path_output_folder,
                    pn.Row(
                        pn.Column(
                            self.predefined_exp_to_cond_title,
                            self.predefined_exp_to_cond
                        ),
                        pn.Card(
                            self.df_exp_to_cond,
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
            height=610,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        self.activate_after_analysis_file_upload = pn.depends(
            self.path_analysis_file.param.value,
            watch=True
        )(self.activate_after_analysis_file_upload)

        self.update_exp_to_cond_df = pn.depends(
            self.predefined_exp_to_cond.param.value,
            watch=True
        )(self.update_exp_to_cond_df)

        self.add_conditions_for_assignment = pn.depends(
            self.df_exp_to_cond.param.value,
            watch=True
        )(self.add_conditions_for_assignment)

        self.run_pipeline = pn.depends(
            self.run_pipeline_button.param.clicks,
            watch=True
        )(self.run_pipeline)

        return LAYOUT


    def natural_sort(self, l):
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
        return sorted(l, key = alphanum_key)


    def extract_sample_names(self):
        with open(self.path_analysis_file.value, 'r') as f:
            all_columns = f.readline().split('\t')
            sample_names = [col for col in all_columns if ('Intensity' in col) and (col != 'Intensity')]
        self.df_exp_to_cond.value = pd.DataFrame(
            data={'sample': self.natural_sort(sample_names), 'condition': str()}
        )


    def add_conditions_for_assignment(self, *args):
        unique_condit = self.df_exp_to_cond.value.condition.unique()
        comb_condit = ['_vs_'.join(comb) for comb in permutations(unique_condit, 2)]
        self.assign_cond_pairs.options = comb_condit


    def set_default_output_folder(self):
        if not self.path_output_folder.value:
            self.path_output_folder.value = os.path.dirname(self.path_analysis_file.value)


    def activate_after_analysis_file_upload(self, *args):
        self.set_default_output_folder()
        self.extract_sample_names()


    def update_exp_to_cond_df(self, *args):
        self.df_exp_to_cond.value = pd.read_csv(
            StringIO(str(self.predefined_exp_to_cond.value, "utf-8")),
            sep='\t'
        )


    def run_pipeline(self, *args):
        import alphaquant.diff_analysis_manager as diffmgr

        self.run_pipeline_progress.active = True

        diffmgr.run_pipeline(
            peptides_tsv=self.path_analysis_file.value,
            samplemap_tsv= StringIO(
                str(self.predefined_exp_to_cond.value, "utf-8")
            ),
            outdir=self.path_output_folder.value,
            pepheader= "peptide",
            protheader="protein"
        )

        self.run_pipeline_progress.active = False
        self.visualize_data_button.clicks += 1
        self.updated.value += 1
