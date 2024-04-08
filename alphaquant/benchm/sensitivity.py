import numpy as np 
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


class RatioClassificationTableGenerator():
    def __init__(self, merged_results_table, decoy_organism,fdr_threshold = 0.05,method_suffixes=["_alphaquant", "_spectronaut"]):
        """This takes in a table that has fdr scored differential expression results from different methods and generates an 
        output table with relevant information for comparing the results of the different methods.

        Example input table columns: "protein" "fdr_alphaquant" "fdr_spectronaut"	"organism_alphaquant"	 "organism_spectronaut" (the organisms need to be specified per method in order to determine the max num allowed decoy hits)
        """
        
        self._merged_results_table = merged_results_table

        self._fdr_column_name = "fdr"
        self._organism_column_name = "organism"
        self._decoy_organism = decoy_organism
        self._method_suffixes = method_suffixes
        self._fdr_threshold = fdr_threshold

        self._per_suffix_results_series = {}

        self.per_species_results_df = pd.DataFrame()

        self._define_per_suffix_results_series()
        self._merge_per_species_results_series_to_df()

    def _define_per_suffix_results_series(self):
        """
        Defines a series of results per method_suffix (method), containing counts of significant results per organism.
        """
        for method_suffix in self._method_suffixes:
            method_results_df = self._merged_results_table[[self._fdr_column_name + method_suffix, self._organism_column_name + method_suffix]].copy()
            hits_per_organism_dict = self._get_hits_per_organism(method_results_df, method_suffix)
            max_hits_per_organism_dict = self._get_max_possible_hits(method_results_df, method_suffix)
            max_allowed_decoy_hits = self._get_max_allowed_decoy_hits(method_results_df, method_suffix)
            
            self._per_suffix_results_series[f"hits{method_suffix}"] = hits_per_organism_dict
            self._per_suffix_results_series[f"max_hits{method_suffix}"] = max_hits_per_organism_dict
            self._per_suffix_results_series[f"allowed_decoy_hits{method_suffix}"] = max_allowed_decoy_hits


    def _get_hits_per_organism(self, suffix_results_df, method_suffix):
        suffix_results_df_significant = self._get_significant_hits(suffix_results_df, method_suffix)
        hits_per_organism = suffix_results_df_significant[self._organism_column_name + method_suffix].value_counts().to_dict()
        return hits_per_organism
    
    def _get_max_possible_hits(self,  suffix_results_df, method_suffix):
        num_entrys_dict = suffix_results_df[self._organism_column_name + method_suffix].value_counts().to_dict()
        num_entrys_dict[self._decoy_organism] = 0
        return num_entrys_dict
    
    def _get_significant_hits(self, suffix_results_df, method_suffix):
        is_significant = suffix_results_df[self._fdr_column_name + method_suffix] < self._fdr_threshold
        return suffix_results_df[is_significant]

    def _get_max_allowed_decoy_hits(self, suffix_results_df, method_suffix):
        suffix_results_df_significant = self._get_significant_hits(suffix_results_df, method_suffix)
        max_allowed_decoy_hits = {} #set non-decoy organisms to nan
        all_organisms = suffix_results_df[self._organism_column_name + method_suffix].unique()
        non_decoy_organisms = [x for x in all_organisms if x != self._decoy_organism]
        max_num_FP = self._fdr_threshold/(1-self._fdr_threshold) * len(suffix_results_df_significant.index)
        max_allowed_decoy_hits[self._decoy_organism] = int(max_num_FP)
        for organism in non_decoy_organisms:
            max_allowed_decoy_hits[organism] = np.nan
        return max_allowed_decoy_hits

    def _merge_per_species_results_series_to_df(self):

        for method_suffix, results_series in self._per_suffix_results_series.items():
            method_df = pd.DataFrame.from_dict(results_series, orient='index', columns=[method_suffix])
            if self.per_species_results_df.empty:
                self.per_species_results_df = method_df
            else:
                self.per_species_results_df = self.per_species_results_df.merge(method_df, left_index=True, right_index=True, how='outer')
        self.per_species_results_df = self.per_species_results_df.reset_index().rename(columns={"index" : "organism"})


def plot_sighits_barplot(df, suffixes, decoy_organism):
    fig, ax = plt.subplots(figsize=(15, 6))

    # Choose a seaborn palette
    palette = sns.color_palette('deep', len(suffixes)) 

    organisms = df['organism']
    n_organisms = len(organisms)
    bar_width = 0.35
    opacity = 0.8
    index = np.arange(n_organisms)

    for i, suffix in enumerate(suffixes):
        hits_col = f'hits{suffix}'
        max_hits_col = f'max_hits{suffix}'

        # Basic color for this suffix
        base_color = palette[i]

        # Slightly modify base color for max hits (lighter)
        max_hits_color = sns.light_palette(base_color, n_colors=3)[1]

        # Plot max hits bars (background)
        ax.bar(index + i * bar_width, df[max_hits_col], bar_width, alpha=0.4, color=max_hits_color, label=max_hits_col)

        # Overlay actual hits bars with slightly darker color
        hits_color = sns.dark_palette(base_color, n_colors=3)[2]
        ax.bar(index + i * bar_width, df[hits_col], bar_width, alpha=opacity, color=hits_color, label=hits_col)

    # Add horizontal lines for allowed decoy hits, if applicable
    if decoy_organism in organisms.values:
        for j, suffix in enumerate(suffixes):
            decoy_value = df.loc[df['organism'] == decoy_organism, f'allowed_decoy_hits{suffix}'].values[0]
            if not np.isnan(decoy_value):
                # Use the base color for the line to match the suffix's color scheme
                line_color = palette[j]
                ax.axhline(y=decoy_value, color=line_color, linestyle='--', label=f'allowed_decoy_hits{suffix}')

    ax.set_xlabel('Organism')
    ax.set_ylabel('Hits')
    ax.set_title('Comparison of Hits by Organism')
    ax.set_xticks(index + bar_width / len(suffixes))
    ax.set_xticklabels(organisms, rotation=45, ha="right")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    fig.tight_layout()
    plt.show()


def get_tp_fp_from_count_df(per_species_results_df, organism_fp, suffix):
    tp_hits = per_species_results_df[per_species_results_df["organism"] != organism_fp][f"hits{suffix}"].sum()
    fp_hits = per_species_results_df[per_species_results_df["organism"] == organism_fp][f"hits{suffix}"].sum()

    return tp_hits, fp_hits
