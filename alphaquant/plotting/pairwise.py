import alphaquant.norm.normalization as aqnorm
import matplotlib.pyplot as plt
import alphaquant.diffquant.diffutils as aq_diff_utils
import numpy as np
import pandas as pd




def plot_normalization_overview(normed_df, samplemap_df):
    normed_df, sample2cond = aq_diff_utils.prepare_loaded_tables(normed_df, samplemap_df)
    sample2cond = dict(zip(samplemap_df["sample"], samplemap_df["condition"]))
    conditions = list(set([sample2cond.get(x) for x in normed_df.columns]))
    conditions = [x for x in conditions if x is not None]
    df_c1 = normed_df[[x for x in normed_df.columns if sample2cond.get(x) == conditions[0]]]
    df_c2 = normed_df[[x for x in normed_df.columns if sample2cond.get(x) == conditions[1]]]

    plot_betweencond_fcs(df_c1, df_c2, True)
    plot_betweencond_fcs(df_c1, df_c2, False)


def plot_withincond_normalization(df_c1, df_c2):
    print("without missingvals (if applicable)")
    plot_betweencond_fcs(aqnorm.drop_nas_if_possible(df_c1), aqnorm.drop_nas_if_possible(df_c2), True)
    print("complete dataset")
    plot_betweencond_fcs(df_c1, df_c2, True)


def plot_betweencond_fcs(df_c1_normed, df_c2_normed, merge_samples=True, cumulative=False):
    """takes normalized intensity dataframes of each condition and plots the distribution of direct peptide fold changes between conditions"""

    if merge_samples:  # samples can be merged to median intensity
        df_c1_normed = df_c1_normed.median(axis=1, skipna=True).to_frame()
        df_c2_normed = df_c2_normed.median(axis=1, skipna=True).to_frame()

    both_idx = df_c1_normed.index.intersection(df_c2_normed.index)
    df1 = df_c1_normed.loc[both_idx]
    df2 = df_c2_normed.loc[both_idx]

    fig, axes = plt.subplots()  # Create a new figure and axes

    for col1 in df1.columns:
        for col2 in df2.columns:
            diff_fcs = df1[col1].to_numpy() - df2[col2].to_numpy()  # calculate fold changes by subtracting log2 intensities of both conditions

            axes.axvline(0, color='red', linestyle="dashed")  # the data is normalized around 0, draw in helper line
            cutoff = max(abs(np.nanquantile(diff_fcs, 0.025)), abs(np.nanquantile(diff_fcs, 0.975)))  # determine 2.5% - 97.5% interval, i.e. remove extremes

            axes.hist(diff_fcs, 80, density=True, histtype='step', range=(-cutoff, cutoff), cumulative=cumulative)  # set the cutoffs to focus the visualization

    axes.set_xlabel("log2(fc)")

    plt.show()
    return fig, axes


def plot_sample_vs_median_fcs(df_c1_normed, df_c2_normed, cumulative=False):
    """Plots the distribution of fold changes between each sample and the median across all samples."""

    # Calculate the median across all samples from both conditions
    combined_median = pd.concat([df_c1_normed, df_c2_normed], axis=1).median(axis=1, skipna=True)

    fig, axes = plt.subplots()  # Create a new figure and axes

    # Compare each sample against the combined median and plot
    for df in [df_c1_normed, df_c2_normed]:
        for col in df.columns:
            diff_fcs = df[col].subtract(combined_median)

            axes.axvline(0, color='red', linestyle="dashed")  # helper line at 0
            cutoff = max(abs(np.nanquantile(diff_fcs, 0.025)), abs(np.nanquantile(diff_fcs, 0.975)))  # determine 2.5% - 97.5% interval

            axes.hist(diff_fcs, 80, density=True, histtype='step', range=(-cutoff, cutoff), cumulative=cumulative)

    axes.set_xlabel("log2(fc)")
    plt.show()
    return fig, axes



def volcano_plot(result_df, fc_header="log2fc", fdr_header="fdr", significance_cutoff=0.05, log2fc_cutoff=0.5, ybound=None, xbound=None):
    result_df[fdr_header] = result_df[fdr_header].replace(0, np.min(result_df[fdr_header].replace(0, 1.0)))
    fdrs = result_df[fdr_header].to_numpy()
    fcs = result_df[fc_header].to_numpy()
    sighits_down = sum((fdrs < significance_cutoff) & (fcs <= -log2fc_cutoff))
    sighits_up = sum((fdrs < significance_cutoff) & (fcs >= log2fc_cutoff))

    fig, ax = plt.subplots()
    ax.set_title(f"{sighits_up} up, {sighits_down} down of {len(fcs)}")

    # Calculate dynamic alpha value

    alpha = max(0.1, min(0.7, 0.7 - 0.6 * (len(fdrs) / 1000)))

    ax.scatter(result_df[fc_header], -np.log10(result_df[fdr_header]), s=10, c='grey', alpha=alpha)

    ax.set_xlabel('log2 FC', fontsize=14)
    ax.set_ylabel('-log10 FDR', fontsize=14)

    if ybound is None:
        ax.set_ylim(0, max(-np.log10(result_df[fdr_header])) + 0.5)
    else:
        ax.set_ylim(ybound)

    if significance_cutoff > 0:
        ax.axhline(y=-np.log10(significance_cutoff), color='g', linestyle='-')
    if log2fc_cutoff > 0:
        ax.axvline(x=log2fc_cutoff, color='g', linestyle='-')
        ax.axvline(x=-log2fc_cutoff, color='g', linestyle='-')

    maxfc = max(abs(result_df[fc_header])) + 0.5
    if xbound is None:
        ax.set_xlim(-maxfc, maxfc)
    else:
        ax.set_xlim(xbound)

    plt.show()
    return fig, ax
