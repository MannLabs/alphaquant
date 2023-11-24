import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_proteoform_intensity_profiles(proteoform_df : pd.DataFrame):
    """creates a panel of proteoform intensity profiles for each protein in the dataframe

    Args:
        proteoform_df (pd.DataFrame): proteoform dataframe loaded from the AlphaQuant output file `medianref_proteoforms.tsv`

    Returns:
        fig, axes: matplotlib figure and axes objects
    """
    fixed_columns = ['proteoform_id', 'peptides', 'number_of_peptides', 'protein', 'corr_to_ref', 'is_reference']
    conditions = proteoform_df.columns.difference(fixed_columns)

    # Group by protein
    grouped = proteoform_df.groupby('protein')
    n_plots = len(grouped)

    # Create a figure with a grid of subplots
    n_cols = int(np.ceil(np.sqrt(n_plots)))
    n_rows = int(np.ceil(n_plots / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))

    # Flatten the axes array for easy indexing
    axes = axes.flatten()

    # Loop through each group and create a plot
    for ax, (protein, sub_df) in zip(axes, grouped):
        for idx, row in sub_df.iterrows():
            style = '-' if row['is_reference'] else '--'
            ax.plot(conditions, row[conditions], style, label=row['proteoform_id'])

        ax.set_title(f'{protein}')
        ax.set_xlabel('Tissue')
        ax.set_ylabel('Expression Level')
        ax.legend(title='Proteoform ID')

    # Turn off any unused subplots
    for i in range(n_plots, len(axes)):
        axes[i].axis('off')
    
    fig.tight_layout()
    plt.show()
    return fig, axes