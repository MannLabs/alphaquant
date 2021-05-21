import os
import re
from io import StringIO
import itertools
import pandas as pd
import numpy as np
from scipy import stats

# alphaquant important
import alphaquant.diffquant_utils as aqutils
import alphaquant.diff_analysis_manager as diffmgr
import alphaquant.visualizations as aqplot

# visualization
import panel as pn
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import seaborn as sns
from holoviews.streams import Selection1D


class BaseWidget(object):

    def __init__(self, name):
        self.name = name
        self.__update_event = pn.widgets.IntInput(value=0)
        self.depends = pn.depends(self.__update_event.param.value)
        self.active_depends = pn.depends(
            self.__update_event.param.value,
            watch=True
        )

    def trigger_dependancy(self):
        self.__update_event.value += 1


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

    def __init__(
        self,
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
            link_url='https://www.biochem.mpg.de/mann',
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

    def __init__(
        self,
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


class RunPipeline(BaseWidget):

    def __init__(self):
        super().__init__(name="Data")
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

    def create(self):
        self.layout = pn.Card(
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
        self.path_analysis_file.param.watch(
            self.activate_after_analysis_file_upload,
            'value'
        )
        self.samplemap.param.watch(
            self.update_samplemap,
            'value'
        )
        self.samplemap_table.param.watch(
            self.add_conditions_for_assignment,
            'value'
        )
        self.run_pipeline_button.param.watch(
            self.run_pipeline,
            'clicks'
        )
        self.visualize_data_button.param.watch(
            self.visualize_data,
            'clicks'
        )
        return self.layout

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
        comb_condit = ['_vs_'.join(comb) for comb in itertools.permutations(unique_condit, 2)]
        self.assign_cond_pairs.options = comb_condit

    def run_pipeline(self, *args):
        self.run_pipeline_progress.active = True

        data_processed, samplemap_df_processed = aqutils.prepare_loaded_tables(
            self.data,
            self.samplemap_table.value
        )

        if self.assign_cond_pairs.value:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.assign_cond_pairs.value]
        else:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.assign_cond_pairs.options]

        diffmgr.run_pipeline(
            unnormed_df=data_processed,
            labelmap_df=samplemap_df_processed,
            results_dir=self.path_output_folder.value,
            condpair_combinations=cond_combinations
        )

        self.trigger_dependancy()
        self.run_pipeline_progress.active = False

    def visualize_data(self, *args):
        self.trigger_dependancy()


class Tabs(object):

    def __init__(self, pipeline):
        self.layout = None
        self.pipeline = pipeline

    def create(
        self,
        # tab_list=None
    ):
        # self.tabs = tab_list
        return self.pipeline.depends(self.create_layout)

    def create_layout(self, *args):
        if self.pipeline.path_output_folder.value is not None and self.pipeline.visualize_data_button.clicks != 0:
            self.pipeline.layout.collapsed = True
            self.layout = pn.Tabs(
                tabs_location='above',
                margin=(30, 10, 5, 8),
                sizing_mode='stretch_width',
            )
            # self.layout += self.tabs
            self.layout += [
                ('Multiple Comparison', MultipleComparison(
                    self.pipeline.path_output_folder.value
                    ).create()
                ),
                ('Single Comparison', SingleComparison(
                    self.pipeline.path_output_folder.value,
                    self.pipeline.samplemap_table.value,
                    ).create()
                ),
            ]
            self.active = 0
            return self.layout


class MultipleComparison(object):

    def __init__(self, output_folder):
        self.condpairs_to_compare = pn.widgets.CrossSelector(
            width=870,
            height=300,
            align='center',
            margin=(15, 15, 15, 15)
        )
        self.output_folder = output_folder
        self.layout = None
        self.heatmap = None

    def create(self):
        self.extract_conditions_from_folder()
        self.condpairs_to_compare.param.watch(
            self.return_clustered_heatmap,
            'value'
        )
        self.layout = pn.Column(
            self.condpairs_to_compare,
            self.heatmap
        )
        return self.layout

    def extract_conditions_from_folder(self):
        self.condpairs_to_compare.options = [f.replace(".results.tsv", "").replace('VS', 'vs') for f in os.listdir(self.output_folder) if re.match(r'.*results.tsv', f)]

    def return_clustered_heatmap(self, *args):
        if self.condpairs_to_compare.value:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.condpairs_to_compare.value]
        else:
            cond_combinations = [tuple(pair.split('_vs_')) for pair in self.condpairs_to_compare.options]

        overview_dataframe = aqplot.get_sample_overview_dataframe(
            results_folder=self.output_folder,
            condpairs_to_compare=cond_combinations
        )

        clustered_dataframe = aqplot.get_clustered_dataframe(overview_dataframe)

        self.layout[1] = pn.Pane(
            self.plot_heatmap(
                clustered_dataframe.T,
                title='Significant proteins heatmap',
                colormap='RdBu'
            )
        )

    def plot_heatmap(self, df, title, colormap):
        """
        This function plots a simple Heatmap.
        """
        print('inside plotting function')
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


class SingleComparison(object):

    def __init__(self, output_folder, sample_to_cond):
        self.condpairs_selector = pn.widgets.Select(
            name='Select a pair of conditions:',
            options=['No conditions'],
            width=870,
            align='center',
            margin=(15, 15, 15, 15)
        )
        self.protein = pn.widgets.AutocompleteInput(
            name='Select protein:',
            placeholder="Type the first letters of the protein's name",
            min_characters=2,
            disabled=True,
            margin=(0, 40, 40, 10)
        )
        self.output_folder = output_folder
        self.sample_to_cond = sample_to_cond
        self.layout = None
        self.volcano_plot = None

    def create(self):
        self.extract_conditions_from_folder()
        self.condpairs_selector.param.watch(
            self.return_volcano_plot,
            'value'
        )
        self.protein.param.watch(
            self.visualize_after_protein_selection,
            'value'
        )
        self.layout = pn.Column(
            self.condpairs_selector,
            None,
            self.volcano_plot,
            pn.layout.Divider(sizing_mode='stretch_width'),
            pn.Row(
                self.protein,
                None,
                None
            ),
            sizing_mode='stretch_width',
            align='center'
        )

        return self.layout

    def extract_conditions_from_folder(self):
        self.condpairs_selector.options.extend([f.replace(".results.tsv", "").replace('VS', 'vs') for f in os.listdir(self.output_folder) if re.match(r'.*results.tsv', f)])

    def return_volcano_plot(self, *args):
        if self.condpairs_selector.value != 'No conditions':
            self.cond1, self.cond2 = self.condpairs_selector.value.split('_vs_')

            self.result_df = aqplot.get_diffresult_dataframe(
                self.cond1,
                self.cond2,
                results_folder=self.output_folder
            )
            self.normalized_intensity_df = aqplot.get_normed_peptides_dataframe(
                self.cond1,
                self.cond2,
                results_folder=self.output_folder
            )
            self.sample2cond = dict(zip(self.sample_to_cond["sample"], self.sample_to_cond["condition"]))

            c1_normed = aqplot.subset_normed_peptides_df_to_condition(
                self.cond1,
                self.sample_to_cond,
                self.normalized_intensity_df
            )
            c2_normed = aqplot.subset_normed_peptides_df_to_condition(
                self.cond2,
                self.sample_to_cond,
                self.normalized_intensity_df
            )

            self.protein.options = self.result_df.protein.values.tolist()
            self.protein.disabled = False

            self.volcano_plot = pn.pane.Plotly(
                self.plot_volcano(
                    self.result_df
                    ),
                align='center',
                config={'responsive': True}
            )

            self.volcano_plot.param.watch(
                self.update_selected_protein,
                'click_data'
            )

            self.layout[1] = pn.Card(
                pn.Pane(
                    self.plot_withincond_fcs(
                        self.normalized_intensity_df,
                        'Log2 fold change distribution between all samples'
                    ),
                    align='center'

                ),
                pn.Pane(
                    self.plot_withincond_fcs(
                        c1_normed,
                        f'Log2 fold change distribution between samples in condition {self.cond1}'
                    ),
                    align='center'
                ),
                pn.Pane(
                    self.plot_withincond_fcs(
                        c2_normed,
                        f'Log2 fold change distribution between samples in condition {self.cond2}'
                    ),
                    align='center'
                ),
                pn.Pane(
                    self.plot_betweencond_fcs(
                        c1_normed,
                        c2_normed
                    ),
                    height=200,
                    width=400,
                    align='center'
                ),
                header='Log2 fold change distribution plots',
                collapsed=True,
                align='center',
                sizing_mode='stretch_width'
            )
            self.layout[2] = self.volcano_plot

        else:
            self.layout[1] = None
            self.layout[2] = None

    def update_selected_protein(self, event):
        self.protein.value = str(event[2].click_data['points'][0]['text'])

    def visualize_after_protein_selection(self, *args):
        result_df_protein_index = self.result_df.set_index("protein")

        melted_df, protein_df = aqplot.get_melted_protein_ion_intensity_table(
            self.protein.value,
            result_df_protein_index,
            self.normalized_intensity_df,
            self.sample2cond
        )

        #get fold change data table for protein
        fc_df = aqplot.get_betweencond_fcs_table(
            melted_df,
            self.cond1,
            self.cond2
        )

        self.layout[4][1] = pn.Pane(
            self.beeswarm_ion_plot(melted_df, protein_df),
            margin=(10, 0, 50, 0)
        )
        self.layout[4][2] = pn.Pane(
            self.foldchange_ion_plot(fc_df, protein_df),
            margin=(10, 0, 50, 0)
        )



    def plot_volcano(
        self,
        result_df,
        fc_header = "log2fc",
        fdr_header = "fdr",
        significance_cutoff = 0.05,
        log2fc_cutoff = 0.5,
        ybound = None,
        xbound = None,
        color='darkgrey',
        marker_size=5,
        name=None,
        opacity=0.9,
        marker_symbol='circle'
    ):
        result_df[fdr_header] = result_df[fdr_header].replace(0, np.min(result_df[fdr_header].replace(0, 1.0)))
        sighits_down = sum((result_df[fdr_header]<significance_cutoff) & (result_df[fc_header] <= -log2fc_cutoff))
        sighits_up = sum((result_df[fdr_header]<significance_cutoff) & (result_df[fc_header] >= log2fc_cutoff))
        result_df_significant = result_df[
            ((result_df[fdr_header] < significance_cutoff) & (result_df[fc_header] <= -log2fc_cutoff)) |
            ((result_df[fdr_header] < significance_cutoff) & (result_df[fc_header] >= log2fc_cutoff))
        ]
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                name='',
                x=result_df[fc_header],
                y=result_df['-log10fdr'],
                mode='markers',
                text=result_df['protein'],
                marker=dict(
                    size=marker_size,
                    symbol=marker_symbol,
                    color=color,
                    opacity=opacity,
                    line=dict(
                        width=1,
                        color='#202020'
                    ),
                    showscale=False
                ),
                hovertemplate =
                    '<b>protein:</b> %{text}'
                    '<br><b>log2fc</b>: %{x:.2f}'+
                    '<br><b>-log10fdr</b>: %{y:.2f}<br>',
            )
        )
        fig.add_trace(
            go.Scatter(
                name='',
                x=result_df_significant[fc_header],
                y=result_df_significant['-log10fdr'],
                mode='markers',
                text=result_df_significant['protein'],
                marker=dict(
                    size=marker_size,
                    symbol=marker_symbol,
                    color='darkgreen',
                    opacity=opacity,
                    line=dict(
                        width=1,
                        color='#202020'
                    ),
                    showscale=False
                ),
                hovertemplate =
                    '<b>protein:</b> %{text}'
                    '<br><b>log2fc</b>: %{x:.2f}'+
                    '<br><b>-log10fdr</b>: %{y:.2f}<br>',
            )
        )
        fig.add_hline(
            y=-np.log10(significance_cutoff),
            line_width=1,
            line_dash="dash",
            line_color="green"
        )
        fig.add_vline(
            x=log2fc_cutoff,
            line_width=1,
            line_dash="dash",
            line_color="green"
        )
        fig.add_vline(
            x=-log2fc_cutoff,
            line_width=1,
            line_dash="dash",
            line_color="green"
        )

        maxfc = max(abs(result_df[fc_header])) + 0.5
        fig.update_layout(
            height=500,
            width=870,
            template='plotly_white',
            title=dict(
                text=f"{sighits_down} down, {sighits_up} up of {len(result_df)}",
                y=0.92,
                x=0.5,
                xanchor='center',
                yanchor='middle',
            ),
            titlefont=dict(size=14, color='black', family='Arial, sans-serif'),
            hovermode='closest',
            xaxis=dict(
                range=[-maxfc,maxfc],
                title='log2 fold change',
                titlefont=dict(size=14, color='black', family='Arial, sans-serif'),
            ),
            yaxis=dict(
                title='-log10 FDR',
                titlefont=dict(size=14, color='black', family='Arial, sans-serif'),
                range=[-0.1, max(-np.log10(result_df[fdr_header])) + 0.5],
            ),
            showlegend=False,
        )

        return fig


    def plot_withincond_fcs(
        self,
        normed_intensity_df,
        title,
        cut_extremes = True
    ):
        """takes a normalized intensity dataframe and plots the fold change distribution between all samples. Column = sample, row = ion"""
        fig = go.Figure()
        samplecombs = list(itertools.combinations(normed_intensity_df.columns, 2))
        cutoff_common = 0
        for spair in samplecombs:#compare all pairs of samples
            s1 = spair[0]
            s2 = spair[1]
            diff_fcs = normed_intensity_df[s1].to_numpy() - normed_intensity_df[s2].to_numpy() #calculate fold changes by subtracting log2 intensities of both samples

            if cut_extremes:
                cutoff_cur = max(abs(np.nanquantile(diff_fcs,0.025)), abs(np.nanquantile(diff_fcs, 0.975))) #determine 2.5% - 97.5% interval, i.e. remove extremes
                if cutoff_cur and cutoff_cur > cutoff_common:
                    cutoff_common = cutoff_cur

            n, bins, _ = plt.hist(diff_fcs, 100, density=True, histtype='step')
            fig.add_trace(
                go.Scatter(
                    x=bins,
                    y=n,
                    mode='lines',
                    line=dict(
                        shape='hvh'
                    )
                )
            )
        fig.update_layout(
            xaxis = dict(
                title="log2 peptide fcs",
                range=(-cutoff_common, cutoff_common)
            ),
            barmode='stack',
            template='plotly_white',
            showlegend=False,
            title = dict(
                    text=title,
                    font=dict(size=16, color='black', family='Arial, sans-serif'),
                    y=0.91,
                    x=0.5,
                    xanchor='center',
                    yanchor='middle',
                ),
            height=400,
            width=870,
        )
        return fig

    def beeswarm_ion_plot(self, df_melted, diffresults_protein, saveloc = None):
        """takes pre-formatted long-format dataframe which contains all ion intensities for a given protein.
          Columns are "ion", "intensity", "condition". Also takes results of the protein differential analysis as a series
          to annotate the plot"""

        fig = plt.figure()

        #get annotations from diffresults
        fdr = float(diffresults_protein.at["fdr"])
        protein = diffresults_protein.name

        #define greyscale color palette for the two conditions
        pal2 = [(0.94, 0.94, 0.94),(1.0, 1.0, 1.0)]

        #searborn standard functions
        ax = sns.boxplot(x="ion", y="intensity", hue="condition", data=df_melted, palette=pal2)
        ax = sns.stripplot(x="ion", y="intensity", hue="condition", data=df_melted, palette="Set2", dodge=True)#size = 10/len(protein_df.index)

        #annotate and format
        handles, labels = ax.get_legend_handles_labels()

        l = plt.legend(handles[2:4], labels[2:4])

        plt.xticks()
        if "gene" in diffresults_protein.index:
            gene = diffresults_line.at["gene"]
            plt.title(f"{gene} ({protein}) FDR: {fdr:.1e}")
        else:
            plt.title(f"{protein} FDR: {fdr:.1e}")

        return fig


    def foldchange_ion_plot(self, df_melted, diffresults_protein, saveloc = None):
        """takes pre-formatted long-format dataframe which contains all between condition fold changes. All ions of a given protein
        are visualized, the columns are "ion" and "log2fc".  Also takes results of the protein differential analysis as a series
          to annotate the plot"""

        fig = plt.figure()

        #get annotations from diffresults
        fdr = float(diffresults_protein.at["fdr"])
        protein = diffresults_protein.name

        #define greyscale color palette
        pal2 = [(0.94, 0.94, 0.94),(1.0, 1.0, 1.0)]

        #plot with seaborn standard functions
        ax = sns.boxplot(x="ion", y="log2fc", data=df_melted, color = "white")
        ax = sns.stripplot(x="ion", y="log2fc", data=df_melted, color = "grey")

        #annotate and format
        handles, labels = ax.get_legend_handles_labels()
        ax.axhline(y = 0, color='black', linewidth=2, alpha=.7, linestyle = "dashed")
        l = plt.legend(handles[2:4], labels[2:4])
        plt.xticks()
        if "gene" in diffresults_protein.index:
            gene = diffresults_line.at["gene"]
            plt.title(f"{gene} ({protein}) FDR: {fdr:.1e}")
        else:
            plt.title(f"{protein} FDR: {fdr:.1e}")

        return fig


    def plot_betweencond_fcs(self, df_c1_normed, df_c2_normed, merge_samples = True):
        """takes normalized intensity dataframes of each condition and plots the distribution of direct peptide fold changes between conditions"""
        fig = plt.figure()
        if merge_samples: #samples can be merged to median intensity
            df_c1_normed = df_c1_normed.median(axis = 1, skipna = True).to_frame()
            df_c2_normed = df_c2_normed.median(axis = 1, skipna = True).to_frame()

        both_idx = df_c1_normed.index.intersection(df_c2_normed.index)
        df1 = df_c1_normed.loc[both_idx]
        df2 = df_c2_normed.loc[both_idx]

        for col1 in df1.columns:
            for col2 in df2.columns:
                diff_fcs = df1[col1].to_numpy() - df2[col2].to_numpy() #calculate fold changes by subtracting log2 intensities of both conditions

                plt.axvline(0, color = 'red', linestyle = "dashed") #the data is normalized around 0, draw in helper line
                cutoff = max(abs(np.nanquantile(diff_fcs,0.025)), abs(np.nanquantile(diff_fcs, 0.975))) #determine 2.5% - 97.5% interval, i.e. remove extremes

                plt.hist(diff_fcs,80,density=True, histtype='step', range=(-cutoff,cutoff)) #set the cutoffs to focus the visualization
        plt.xlabel("log2(fc)")

        return fig
