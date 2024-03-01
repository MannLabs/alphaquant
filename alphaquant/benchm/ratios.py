import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

class MixedSpeciesScatterPlotter():
    """
    Plots the LFQ-bench style plots from a standardized input table. The columns of an example table are:
    'protein'	'log2fc_alphaquant'	'intensity_alphaquant' 'organism'	'log2fc_spectronaut' 'intensity_spectronaut'
    """
    def __init__(self, df_combined, method_suffixes, expected_log2fcs, figure_size = [4, 4], point_size = 4, alpha = 0.5):
        self._df_combined = df_combined
        self._method_suffixes = method_suffixes
        self._expected_log2fcs = expected_log2fcs
        self._figure_size = figure_size
        self._point_size = point_size
        self._alpha = alpha
        

        self.fig = None
        self.axes = None

        self._plot_fc_scatter_per_method()

    
    def _plot_fc_scatter_per_method(self):
        self._prepare_fig_and_axes()
        for method_idx in range(len(self._method_suffixes)):
            self._plot_fc_scatter(method_idx)
        self._set_uniform_axis_ranges()
        self._create_unified_legend()
        self.fig.tight_layout()

    def _prepare_fig_and_axes(self):
        num_methods = len(self._method_suffixes)
        self.fig, self.axes = plt.subplots(1, num_methods, figsize=(self._figure_size[0]*num_methods, self._figure_size[1]), squeeze=False)
    
    def _plot_fc_scatter(self, method_idx):
        suffix = self._method_suffixes[method_idx]
        intensity_column = f'intensity{suffix}'
        log2fc_column = f'log2fc{suffix}'
        organism_column = 'organism'
        ax = self.axes[0][method_idx]
        sns.scatterplot(data = self._df_combined, x=intensity_column, y=log2fc_column, hue=organism_column, ax=ax, s = self._point_size, alpha = self._alpha)
        ax.set_xscale('log')
        for expected_log2fc in self._expected_log2fcs:
            ax.axhline(expected_log2fc, color='black')
        ax.set_title(self._get_title_text(suffix,log2fc_column, organism_column))
        ax.set_xlabel("intensity")
        ax.set_ylabel("log2 fold change")
        ax.get_legend().remove()

    def _get_title_text(self, suffix,log2fc_column, organism_column):
        std_devs = self._df_combined.groupby(organism_column)[log2fc_column].std()
    
        title_text = f"{suffix[1:]}:\n"
        for organism, std_dev in std_devs.items():
            title_text += f"{organism}: {std_dev:.2f}\n"
        return title_text
    
    def _set_uniform_axis_ranges(self):
        # Find the overall min and max across all subplots for both axes
        all_x_lims = [ax.get_xlim() for ax in self.axes.flatten()]
        all_y_lims = [ax.get_ylim() for ax in self.axes.flatten()]

        global_x_min = min(lim[0] for lim in all_x_lims)
        global_x_max = max(lim[1] for lim in all_x_lims)
        global_y_min = min(lim[0] for lim in all_y_lims)
        global_y_max = max(lim[1] for lim in all_y_lims)

        # Set the same x and y limits for all plots
        for ax in self.axes.flatten():
            ax.set_xlim(global_x_min, global_x_max)
            ax.set_ylim(global_y_min, global_y_max)
    
    def _create_unified_legend(self):
        # Only create a unified legend if there's more than one category
        if len(set(self._df_combined['organism'])) > 1:
            # Place the legend on the right side of the last subplot
            self.axes[0, -1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Organism')



import plotly.express as px

class MixedSpeciesScatterPlotterInteractive:
    """
    Plots interactive LFQ-bench style plots from a standardized input table using plotly for interactivity.
    Hovering over points will display protein names among other details.
    """
    def __init__(self, df_combined, method_suffixes, expected_log2fcs):
        self._df_combined = df_combined
        self._method_suffixes = method_suffixes
        self._expected_log2fcs = expected_log2fcs
        self.figures = []

        self._plot_fc_scatter_per_method()

    def _plot_fc_scatter_per_method(self):
        for method_idx, suffix in enumerate(self._method_suffixes):
            self._plot_fc_scatter(method_idx, suffix)

    def _plot_fc_scatter(self, method_idx, suffix):
        intensity_column = f'intensity{suffix}'
        log2fc_column = f'log2fc{suffix}'
        organism_column = 'organism'
        protein_column = 'protein'

        fig = px.scatter(self._df_combined, x=intensity_column, y=log2fc_column, color=organism_column,
                          hover_data=[protein_column], log_x=True, title=suffix[1:],
                          labels={"x": "Intensity", "y": "Log2 Fold Change"},
                          opacity=0.5, size_max=60)

        # Add expected_log2fc lines
        for expected_log2fc in self._expected_log2fcs:
            fig.add_hline(y=expected_log2fc, line_dash="dash", line_color="black")

        self.figures.append(fig)
        print("scatter plotted")

    def show_figures(self):
        for fig in self.figures:
            fig.show()



import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

class MixedSpeciesBoxPlotter():
    """
    Plots box plots for log2 fold changes across organisms from a standardized input table. 
    The columns of an example table are:
    'protein', 'log2fc_alphaquant', 'organism', 'log2fc_spectronaut'
    """
    def __init__(self, df_combined, method_suffixes, expected_log2fcs, figure_size = [4, 4]):
        self._df_combined = df_combined
        self._method_suffixes = method_suffixes
        self._expected_log2fcs = expected_log2fcs
        self._figure_size = figure_size


        self.fig = None
        self.axes = None

        self._plot_box_per_method()

    def _plot_box_per_method(self):
        self._prepare_fig_and_axes()
        for method_idx in range(len(self._method_suffixes)):
            self._plot_box(method_idx)
        self._set_uniform_axis_ranges()
        self.fig.tight_layout()

    def _prepare_fig_and_axes(self):
        num_methods = len(self._method_suffixes)
        self.fig, self.axes = plt.subplots(1, num_methods, figsize=(self._figure_size[0] * num_methods, self._figure_size[1]), squeeze=False)

    def _plot_box(self, method_idx):
        suffix = self._method_suffixes[method_idx]
        log2fc_column = f'log2fc{suffix}'
        organism_column = 'organism'
        ax = self.axes[0][method_idx]
        sns.boxplot(data=self._df_combined, x=organism_column, y=log2fc_column, ax=ax)
        for expected_log2fc in self._expected_log2fcs:
            ax.axhline(expected_log2fc, color='black')
        title_text = self._get_title_text(suffix,log2fc_column, organism_column)
        ax.set_title(title_text)
        ax.set_xlabel("organism")
        ax.set_ylabel("log2 Fold Change")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')

    def _get_title_text(self, suffix,log2fc_column, organism_column):
        std_devs = self._df_combined.groupby(organism_column)[log2fc_column].std()
    
        title_text = f"{suffix[1:]}:\n"
        for organism, std_dev in std_devs.items():
            title_text += f"{organism}: {std_dev:.2f}\n"
        return title_text

    def _set_uniform_axis_ranges(self):
        # Collect all y-values from all subplots
        all_y_values = []
        for ax in self.axes.flatten():
            for line in ax.get_lines():
                all_y_values.extend(line.get_ydata())
        
        # Calculate the percentiles to exclude severe outliers
        lower_percentile = 0.5  # Adjust the percentile as needed
        upper_percentile = 100-lower_percentile
        global_y_min = np.percentile(all_y_values, lower_percentile)
        global_y_max = np.percentile(all_y_values, upper_percentile)

        # Set the same y limits for all plots, excluding severe outliers
        for ax in self.axes.flatten():
            ax.set_ylim(global_y_min, global_y_max)



class MixedSpeciesRatioComparison():
    """
    Provides the comparison of observed log2fcs and expected log2fcs The columns of an example table are:
    'protein'	'log2fc_alphaquant'	'intensity_alphaquant' 'organism'	'log2fc_spectronaut' 'intensity_spectronaut'
    """
    def __init__(self, df_combined, method_suffixes, organism2expectedfc):
        self._df_combined = df_combined
        self._method_suffixes = method_suffixes
        self._organism2expectedfc = organism2expectedfc

        self._results_dict_list = []
        self.results_df = None

        self._fill_results_dict_list()
        self.results_df = pd.DataFrame(self._results_dict_list)

    
    def _fill_results_dict_list(self):
        for method_idx in range(len(self._method_suffixes)):
            self._collect_fc_distance_metrics_per_organism(method_idx)

    def _collect_fc_distance_metrics_per_organism(self, method_idx):
        method_suffix = self._method_suffixes[method_idx]
        log2fc_column = f'log2fc{method_suffix}'
        organism_column = 'organism'

        df_organism_grouped = self._df_combined.groupby(organism_column)
        for organism, df_organism in df_organism_grouped:
            organism_fcs = df_organism[log2fc_column].to_numpy()
            self._add_fc_distance_metrics_to_results_dict_list(organism_fcs, organism, method_suffix)
    
    def _add_fc_distance_metrics_to_results_dict_list(self, organism_fcs, organism, method_suffix):
        expected_fc = self._organism2expectedfc[organism]
        fc_distance = np.abs(organism_fcs - expected_fc)
        median_fc_distance = np.median(fc_distance)
        upper_10th_percentile = np.percentile(fc_distance, 90)
        lower_10th_percentile = np.percentile(fc_distance, 10)
        self._results_dict_list.append({'method_suffix': method_suffix,'organism': organism, 'median_fc_distance': median_fc_distance, 
                                        'upper_10th_percentile': upper_10th_percentile, 'lower_10th_percentile': lower_10th_percentile})
        
    
    
