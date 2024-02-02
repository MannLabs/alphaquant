import seaborn as sns
import matplotlib.pyplot as plt

class MixedSpeciesFCPlotter():
    """
    Plots the LFQ-bench style plots from a standardized input table. The columns of an example table are:
    'protein'	'log2fc_alphaquant'	'intensity_alphaquant' 'organism'	'log2fc_spectronaut' 'intensity_spectronaut'
    """
    def __init__(self, df_combined, method_suffixes, expected_log2fcs ):
        self._df_combined = df_combined
        self._method_suffixes = method_suffixes
        self._expected_log2fcs = expected_log2fcs

        self.fig = None
        self.axes = None

        self._prepare_fig_and_axes()
        self._plot_fc_scatter_per_method()

    
    def _prepare_fig_and_axes(self):
        num_methods = len(self._method_suffixes)
        self.fig, self.axes = plt.subplots(1, num_methods, figsize=(5*num_methods, 5))
    
    def _plot_fc_scatter_per_method(self):
        for method_idx in range(len(self._method_suffixes)):
            self._plot_fc_scatter(method_idx)
        self.fig.tight_layout()
    
    def _plot_fc_scatter(self, method_idx):
        suffix = self._method_suffixes[method_idx]
        intensity_column = f'intensity{suffix}'
        log2fc_column = f'log2fc{suffix}'
        organism_column = 'organism'
        ax = self.axes[method_idx]
        sns.scatterplot(data = self._df_combined, x=intensity_column, y=log2fc_column, hue=organism_column, ax=ax)
        for expected_log2fc in self._expected_log2fcs:
            ax.axhline(expected_log2fc)
        ax.set_title(suffix[1:])
        ax.set_xlabel("intensity")
        ax.set_ylabel("log2 fold change")
