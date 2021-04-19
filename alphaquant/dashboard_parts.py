import os
from io import StringIO

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
            placeholder='Enter the whole path to the MQ / Spectronaut / DIA-NN output file',
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
            margin=(0,0,0,12)
        )
        self.predefined_exp_to_cond = pn.widgets.FileInput(
            accept='.tsv,.csv',
            margin=(-10,0,5,12)
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
            margin=(-10, 0, 30, 0)
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
                    self.predefined_exp_to_cond_title,
                    self.predefined_exp_to_cond,
                    margin=(10, 30, 10, 10),
                ),
                pn.Spacer(sizing_mode='stretch_width'),
                pn.Column(
                    self.run_pipeline_button,
                    self.run_pipeline_progress,
                    align='center',
                    margin=(100, 40, 0, 0),
                )
            ),
            title='Data Import / Run Pipeline',
            collapsed=False,
            header_background='#eaeaea',
            header_color='#333',
            align='center',
            sizing_mode='stretch_width',
            height=225,
            margin=(5, 8, 10, 8),
            css_classes=['background']
        )

        self.run_pipeline = pn.depends(
            self.run_pipeline_button.param.clicks,
            watch=True
        )(self.run_pipeline)

        return LAYOUT


    def run_pipeline(self, *args):
        import alphaquant.diff_analysis_manager as diffmgr

        self.run_pipeline_progress.active = True

        diffmgr.run_pipeline(
            peptides_tsv=self.path_analysis_file.value,
            samplemap_tsv= StringIO(
                str(self.predefined_exp_to_cond.value, "utf-8")
            ),
            pepheader= "peptide",
            protheader="protein",
            outdir=self.path_output_folder.value
        )

        self.run_pipeline_progress.active = False
        self.updated.value += 1
