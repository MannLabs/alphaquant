import seaborn as sns
import matplotlib.pyplot as plt

class MixedSpeciesScatterPlotter():
    """
    Plots the LFQ-bench style plots from a standardized input table. The columns of an example table are:
    'protein'	'log2fc_alphaquant'	'intensity_alphaquant' 'organism'	'log2fc_spectronaut' 'intensity_spectronaut'
    """
    def __init__(self, df_combined, method_suffixes, expected_log2fcs, figure_size = [4, 4]):
        self._df_combined = df_combined
        self._method_suffixes = method_suffixes
        self._expected_log2fcs = expected_log2fcs
        self._figure_size = figure_size
        

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
        sns.scatterplot(data = self._df_combined, x=intensity_column, y=log2fc_column, hue=organism_column, ax=ax)
        for expected_log2fc in self._expected_log2fcs:
            ax.axhline(expected_log2fc, color='black')
        ax.set_title(suffix[1:])
        ax.set_xlabel("intensity")
        ax.set_ylabel("log2 fold change")
        ax.get_legend().remove()
    
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
        ax.set_title(suffix[1:])
        ax.set_xlabel("organism")
        ax.set_ylabel("log2 Fold Change")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')

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
