
import alphaquant.diffquant.diffutils as utils
import alphaquant.config.variables as aq_variables

# Cell
import pandas as pd
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt


class AlphaPeptColorMap():
    def __init__(self):
        self.colorlist_hex  = ["#3FC5F0", "#42DEE1", "#7BEDC5", "#FFD479", "#16212B"]
        self.colorlist = [matplotlib.colors.to_rgba(x) for x in self.colorlist_hex]

        self.colormap_linear = matplotlib.colors.LinearSegmentedColormap.from_list("alphapept",self.colorlist)
        self.colormap_discrete = matplotlib.colors.LinearSegmentedColormap.from_list("alphapept",self.colorlist, N=5)

class ClusterColorMap():
    def __init__(self):
        self.colorlist_hex = [    "#D32F2F",  # Crimson Red
                    "#FFA000",  # Burnt Orange
                    "#FFEB3B",  # Golden Yellow
                    "#4CAF50",  # Grass Green
                    "#00BCD4",  # Cyan Blue
                    "#303F9F",  # Cobalt Blue
                    "#7B1FA2",  # Deep Purple
                    "#E91E63",  # Rose Pink
                    "#795548",  # Mocha Brown
                    "#607D8B"   # Slate Grey
                    ]
        self.colorlist = [matplotlib.colors.to_rgba(x) for x in self.colorlist_hex]

class AlphaPeptColorMapAdapted():
    def __init__(self):
        self.colorlist_hex  =  [
    "#3FB7E4",  # Vivid Sky Blue (Slightly desaturated)
    "#7BEDC5",  # Medium Aquamarine
    "#EBCB70",  # Mustard (Slightly less bright)
    "#16212B",  # Gunmetal
    "#CA9ECB",  # Soft Lilac (Slightly more saturated)
    "#708090",  # Slate Gray
    "#3391A6",  # Deep Cerulean (Slightly lightened)
    "#AEDDE9",  # Powder Blue (Slightly warmer)
    "#5F9EA0",  # Cadet Blue
    "#E77D7D"   # Light Coral (Slightly desaturated)
]
        self.colorlist = [matplotlib.colors.to_rgba(x) for x in self.colorlist_hex]

class AlphaQuantColorMap():
    def __init__(self):
        self.colorlist_hex  = ["#d8674e",  # Cadmium Red
        "#45a6ce",  # Steel Blue
        "#fdb73b",  # Cadmium Yellow
        "#a6d1f1",  # Baby Blue
        "#b04e8d",  # Tiffany Rose
        "#6e79b9",  # Periwinkle
        "#fcdf3b",  # Goldenrod
        "#50C878",  # Emerald Green
        "#808080",  # Grey instead of Amber
        "#FF7F50",  # Coral
        "#0F52BA",  # Egyptian Blue
        "#9966CC",  # Amethyst
        "#40E0D0"   # Turquoise
        ]
        self.colorlist = [matplotlib.colors.to_rgba(x) for x in self.colorlist_hex]


def rgba_list_to_hex_list(rgba_list):
    hex_list = []
    for rgba in rgba_list:
        # Convert each value to a 0-255 scale, then to hex, and finally concatenate.
        hex_code = '#' + ''.join([f"{int(c*255):02X}" for c in rgba[:3]])
        hex_list.append(hex_code)
    return hex_list

# Cell
import matplotlib.pyplot as plt
import numpy as np

def plot_pvals(result_df):
    pvals = result_df["peptide_pval"].to_list()
    plt.hist(pvals,99,cumulative=True,density=True, histtype='step')
    x = np.linspace(0,1,100)
    plt.plot(x, x)
    plt.show()

# Cell
import matplotlib.pyplot as plt
def scatter_df_columns(merged_df, log_axes = False):
    col = (0.2, 0.4, 0.6, 0.1)
    ref_columns = list(filter(lambda x : "_ref" in x, merged_df.columns.to_list())) #filter the reference columns from the merged df

    for ref in ref_columns:
        compare = ref.replace("_ref", "")
        ax_p = merged_df.plot.scatter(x=ref,y=compare, color = col)
        corr = merged_df[ref].corr(merged_df[compare])
        plt.title(f"{ref} vs. {compare} corr {corr}")
        x = np.linspace(0,merged_df[ref].max(),100)
        plt.plot(x, x)
        if log_axes:
            plt.xscale('log')
            plt.yscale('log')
        plt.show()

# Cell
import anytree
import time
def get_melted_protein_ion_intensity_table(protein, diffresults_df, normed_df, sample2cond, condpair_root_node = None,ion_header = 'quant_id', protein_header = 'protein'):
    t_start = time.time()
    diffresults_line = diffresults_df.loc[protein]
    value_vars = set.intersection(set(normed_df.columns), set(sample2cond.keys()))
    protein_df = normed_df.xs(protein, level = 0)
    df_melted = pd.melt(protein_df.reset_index(), value_vars= value_vars, id_vars=[ion_header], value_name="intensity", var_name="sample")
    df_melted["condition"] = [sample2cond.get(x) for x in df_melted["sample"]]
    t_melted = time.time()
    #if ion clustering has been performed, add cluster information
    if condpair_root_node != None:

        protein_node = anytree.findall_by_attr(condpair_root_node, protein, maxlevel=2)[0]
        ions_sorted = [x.name for x in protein_node.leaves]
        ion2is_included = {x.name : x.cluster==0 for x in protein_node.leaves} #written as dict because identical ion has multiple columns
        ions_in_df = set(df_melted[ion_header]) - set(ions_sorted)
        if len(ions_in_df)>0:
            Exception("Clustered ions and observed ions differ!")

        df_melted = df_melted.set_index(ion_header)
        df_melted = df_melted.loc[ions_sorted]
        df_melted["is_included"] = [ion2is_included.get(x) for x in df_melted.index]
        df_melted = df_melted.reset_index()
    t_annotated = time.time()
    print(f"times melted protein intensities:\n t_melted: {t_melted - t_start} \n t_annotated: {t_annotated - t_melted}")
    return df_melted, diffresults_line

# Cell
import pandas as pd
import os
import numpy as np

def get_diffresult_dataframe(cond1, cond2, results_folder = os.path.join(".", "results")):
    """
    reads the results dataframe for a given condpair
    """
    condpair = utils.get_condpairname([cond1, cond2])
    diffresults = os.path.join(results_folder, f"{condpair}.results.tsv")

    try:
        diffprots = pd.read_csv(diffresults, sep = "\t")
    except:
        print(f"no quantfiles found for {condpair}!")
        return None
    diffprots = diffprots[(diffprots["condition_pair"] == condpair)]

    diffprots["-log10fdr"] = -np.log10(diffprots["fdr"])
    #diffprots = diffprots.set_index("protein")

    return diffprots

# Cell
import pandas as pd
import os
import numpy as np

def get_normed_peptides_dataframe(cond1, cond2, results_folder = os.path.join(".", "results")):
    condpair = utils.get_condpairname([cond1, cond2])
    normed_peptides_tsv = os.path.join(results_folder, f"{condpair}.normed.tsv")
    try:
        normed_peptides = pd.read_csv(normed_peptides_tsv, sep = "\t")
    except:
        print(f"no normed peptides found for {condpair}!")
        return None

    numeric_cols = list(normed_peptides.select_dtypes(include=np.number).columns)
    #available_vals = list(set(samplemap_df["sample"].values).intersection(set(normed_peptides.columns)))
    normed_peptides[numeric_cols] = np.log2(normed_peptides[numeric_cols].replace(0, np.nan))
    normed_peptides = normed_peptides.set_index(["protein", aq_variables.QUANT_ID])
    return normed_peptides

# Cell
import pandas as pd

def initialize_sample2cond(samplemap):
    samplemap_df = pd.read_csv(samplemap, sep = "\t")
    sample2cond = dict(zip(samplemap_df["sample"], samplemap_df["condition"]))
    return samplemap_df, sample2cond

# Cell
import pandas as pd
import numpy as np
import plotly.graph_objects as go

def plot_volcano_plotly(
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
                color=AlphaQuantColorMap().colorlist_hex[1],
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
        line_color=AlphaQuantColorMap().colorlist_hex[1]
    )
    fig.add_vline(
        x=log2fc_cutoff,
        line_width=1,
        line_dash="dash",
        line_color=AlphaQuantColorMap().colorlist_hex[1]
    )
    fig.add_vline(
        x=-log2fc_cutoff,
        line_width=1,
        line_dash="dash",
        line_color=AlphaQuantColorMap().colorlist_hex[1]
    )

    maxfc = max(abs(result_df[fc_header])) + 0.5
    fig.update_layout(
        height=500,
        width=870,
        template='plotly_white',
        title=dict(
            text=f"{sighits_down} down, {sighits_up} up of {len(result_df)}",
            font=dict(size=14, color='black', family='Arial, sans-serif'),
            y=0.92,
            x=0.5,
            xanchor='center',
            yanchor='middle',
        ),
        hovermode='closest',
        xaxis=dict(
            range=[-maxfc,maxfc],
            title=dict(
                text='log2 fold change',
                font=dict(size=14, color='black', family='Arial, sans-serif'),
            )
        ),
        yaxis=dict(
            title=dict(
                text='-log10 FDR',
                font=dict(size=14, color='black', family='Arial, sans-serif'),
            ),
            range=[-0.1, max(-np.log10(result_df[fdr_header])) + 0.5],
        ),
        showlegend=False,
    )

    return fig

# Cell
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns

def scatter_ml_regression(y_test, y_pred, results_dir = None):
    fig, ax = plt.subplots()

    sns.regplot(x = y_test, y = y_pred, scatter_kws=dict(alpha=0.1), ax = ax)
    err = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    ax.set_title(f"MSE: {err:.2f}, R2: {r2:.2f}")

    if results_dir is not None:
        fig.savefig(f"{results_dir}/ml_regression.pdf")
    plt.show()



def plot_perturbed_unperturbed_fcs(fcs_perturbed, fcs_unperturbed, results_dir = None):
    fig, ax = plt.subplots()
    ax.hist(bins = 60, x=fcs_unperturbed, label= 'unperturbed',density=True, histtype='step')
    ax.hist(bins = 60, x= fcs_perturbed,label='perturbed', density=True, histtype='step')
    if results_dir is not None:
        ax.figure.savefig(f"{results_dir}/compare_pertubed_unperturbed.pdf")


# Cell
import seaborn as sns
import numpy as np

def plot_fc_intensity_scatter(result_df, name, ax = None,expected_log2fc = None, tolerance_interval = 0.5, xlim_lower = -1, xlim_upper = 3.5):
    result_df["log_median_intensity"] = np.log10(result_df["median_intensity"])
    if ax == None:
        ax = sns.scatterplot( x="log2fc", y="log_median_intensity", data=result_df, alpha=0.2)
    else:
        sns.scatterplot(x="log2fc", y="log_median_intensity", data=result_df, alpha=0.2, ax = ax)
    if expected_log2fc is not None:
        ax.vlines(expected_log2fc, 3.8, 11)
        ax.vlines(expected_log2fc-tolerance_interval, 3.8, 11, linestyles = 'dotted')
        ax.vlines(expected_log2fc+tolerance_interval, 3.8, 11, linestyles = 'dotted')
        ax.set(xlim = (xlim_lower, expected_log2fc + xlim_upper))
    std = np.std(result_df['log2fc'].values)
    mean = np.mean(result_df['log2fc'].values)
    ax.set_title(f'{name}\n mean {mean:.2f}, std {std:.2f}, nums {len(result_df["log2fc"])}')
    #ax.set(xlim = (-2, 6))
    #ax.set(ylim = (3.8, 11))
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if ax == None:
        plt.show()

# Cell
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_violin_plots_log2fcs(names, dfs, ax = None):
    df_fcs_longformat = get_longformat_df(names, dfs)
    sns.violinplot(x='variable', y='value', data=df_fcs_longformat, scale='width', ax = ax)
    #sns.stripplot(x='variable', y='value', data=df_fcs_longformat, ax = ax, alpha = 0.35, color = 'lightgrey')
    if ax == None:
        plt.show()

def plot_beeswarm_plot_log2fcs(names, dfs, ax = None):
    df_fcs_longformat = get_longformat_df(names, dfs)
    sns.boxplot(x='variable', y='value', data=df_fcs_longformat, ax = ax)
    sns.stripplot(x='variable', y='value', data=df_fcs_longformat, ax = ax, alpha = 0.35)
    if ax == None:
        plt.show()

def get_longformat_df(names, dfs):
    methods = []
    fcs = []
    for idx in range(len(names)):
        df = dfs[idx]
        name = names[idx]
        fcs_local = list(df["log2fc"])
        fcs.extend(fcs_local)
        methods.extend([name for x in range(len(fcs_local))])
    df_fcs_longformat = pd.DataFrame({'variable' : methods, 'value' : fcs})
    return df_fcs_longformat

# Cell
from matplotlib import pyplot as plt
import numpy as np

def plot_feature_importances(coef, names, top_n = np.inf, print_out_name = False, results_dir = None):
    imp,names = filter_sort_top_n(coef, names, top_n)
    fig, ax = plt.subplots()
    ax.set_title('Feature Importances')
    ax.barh(range(len(names)), imp, align='center')
    ax.set_yticks(range(len(names)), names)
    if results_dir is not None:
        fig.savefig(f"{results_dir}/ml_feature_importances.pdf")
    plt.show()


def filter_sort_top_n(imp, names, top_n):
    tuplelist = list(zip(imp, names))
    tuplelist.sort(key = lambda x : abs(x[0]),reverse=True)
    tuplelist = tuplelist[:top_n]
    imp = [x[0] for x in tuplelist]
    names = [x[1] for x in tuplelist]
    return imp, names

# Cell

import matplotlib.pyplot as plt
import numpy as np
from sklearn import metrics
import random
import seaborn as sns
import alphaquant.diffquant.diffutils as aqdiffutils




def plot_predictability_roc_curve( true_falses, ml_scores, reference_scores, ax = None, percentile_cutoff_indication = None):

    if percentile_cutoff_indication is not None:
        ax.axhline(percentile_cutoff_indication, color = 'lightgrey')
    plot_roc_curve(true_falses, ml_scores, "AlphaQuant score", ax)
    plot_roc_curve(true_falses, reference_scores, f"reference score", ax)
    true_falses = random.sample(true_falses, len(true_falses))
    plot_roc_curve(true_falses, ml_scores, f"random score", ax)
    ax.set_title('ROC curve')
    ax.set_xlabel('FPR')
    ax.set_ylabel('TPR')
    ax.legend()


def plot_predictability_precision_recall_curve( true_falses, ml_scores, reference_scores, ax = None, percentile_cutoff_indication = None):

    if percentile_cutoff_indication is not None:
        ax.axvline(percentile_cutoff_indication, color = 'lightgrey')
    plot_precision_recall_curve(true_falses, ml_scores, "AlphaQuant score", ax)
    plot_precision_recall_curve(true_falses, reference_scores, "reference score", ax)
    true_falses = random.sample(true_falses, len(true_falses))
    plot_precision_recall_curve(true_falses, ml_scores, "random score", ax)
    ax.set_title('Precision-recall curve')
    ax.set_xlabel('recall')
    ax.set_ylabel('precision')
    ax.legend()


def plot_outlier_fraction(node_df, reference_df, expected_log2fc, outlier_thresholds, ax = None):
    thresholds = []
    aq_fractions = []
    reference_fractions = []
    for threshold in outlier_thresholds:
        thresholds.append(threshold)
        aq_fractions.append(aqdiffutils.count_fraction_outliers_from_expected_fc(node_df,threshold, expected_log2fc))
        reference_fractions.append(aqdiffutils.count_fraction_outliers_from_expected_fc(reference_df, threshold,expected_log2fc))
    df = pd.DataFrame({'threshold' :thresholds, 'AlphaQuant' : aq_fractions, 'reference' : reference_fractions})
    df_unpiv = df.melt(id_vars = ['threshold'])
    sns.barplot(x = "threshold", y = 'value', hue = 'variable', data = df_unpiv, ax = ax)







def get_true_false_to_ml_scores(nodes, expected_fc, fc_cutoff_bad = 1, fc_cutoff_good = 0.3, reverse = False):
    true_falses = []
    ml_scores = []
    reference_scores = []
    fcs = []

    for node in nodes:
        fc_diff = abs(node.fc - expected_fc)
        if fc_diff>fc_cutoff_bad:
            true_falses.append(False)
            ml_scores.append(1/abs(node.ml_score))
            reference_scores.append(node.default_quality_score)
            fcs.append(node.fc)
        if fc_diff<fc_cutoff_good:
            true_falses.append(True)
            ml_scores.append(1/abs(node.ml_score))
            reference_scores.append(node.default_quality_score)
            fcs.append(node.fc)

    if reverse:
        true_falses = [not x for x in true_falses]
        ml_scores = [1/x for x in ml_scores]

    print(f"num trues{sum(true_falses)}\tnum falses {len(true_falses) - sum(true_falses)}")

    return true_falses, ml_scores, reference_scores, fcs

def plot_true_false_fcs_of_test_set(fcs, true_falses, ax):
    plot_dict = {'fcs': fcs, 'true_false' : true_falses}
    sns.stripplot(data = plot_dict, x = "true_false", y = 'fcs', ax=ax, palette=[sns.color_palette()[3], sns.color_palette()[0]])
    ax.set_ylabel('log2FC')
    ax.set_xlabel('Cathegory for ROC curve')

def plot_fc_dist_of_test_set(fcs, ax):
    ax.hist(fcs, 60, density=True, histtype='step')
    ax.set_xlabel('log2FC')
    median = np.median(fcs)
    plt.axvline(x=median)
    ax.set_title(f'FC distribution of set, median {median}')



def plot_roc_curve(true_falses, scores, name, ax):
    fpr, tpr, _ = metrics.roc_curve(true_falses,  scores)
    if ax is not None:
        ax.plot(fpr,tpr, label = name)
    else:
        plt.plot(fpr,tpr, label = name)
        plt.legend()
        plt.show()

def plot_precision_recall_curve(true_falses, scores, name, ax):
    precision, recall, _ = metrics.precision_recall_curve(true_falses, scores)
    if ax is not None:
        ax.plot(recall, precision, label = name)
    else:
        plt.plot(recall, precision, label = name)
        plt.legend()
        plt.show()



# Cell
import alphaquant.plotting.base_functions as aqviz
import alphaquant.utils.utils as aqutils

import anytree
import matplotlib.pyplot as plt

def compare_fcs_unperturbed_vs_perturbed_and_clustered(results_dir_unperturbed, results_dir_perturbed, results_dir_perturbed_unclustered):
    ctree_unperturbed = aqutils.read_condpair_tree("S1_filtered", "S2_filtered",results_folder=results_dir_unperturbed)
    ctree_perturbed = aqutils.read_condpair_tree("S1_annot", "S2_annot", results_folder = results_dir_perturbed)
    results_df_perturbed_unclustered =  aqviz.get_diffresult_dataframe("S1", "S2", results_dir_perturbed_unclustered)
    fcs_unperturbed = [x.fc for x in ctree_unperturbed.children]
    fcs_perturbed = [x.fc for x in ctree_perturbed.children]
    fcs_perturbed_unclustered = results_df_perturbed_unclustered["log2fc"]
    plt.hist(fcs_unperturbed,label = "unperturbed", cumulative=True, bins=50, histtype='step')
    plt.hist(fcs_perturbed,label = "perturbed", cumulative=True, bins=50, histtype='step')
    plt.hist(fcs_perturbed_unclustered,label = "perturbed_unclustered", cumulative=True, bins=50, histtype='step')
    plt.legend()
    plt.show()


def rgb_to_hex(rgb):
    if len(rgb) == 3:
        return "#{:02x}{:02x}{:02x}".format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
    elif len(rgb) == 4:
        return "#{:02x}{:02x}{:02x}{:02x}".format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255), int(rgb[3] * 255))
    else:
        raise ValueError("RGB input not recognized")


