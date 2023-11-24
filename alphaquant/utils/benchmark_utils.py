import numpy as np
def subset_df_to_n_most_complete_proteins(proteome_df_aq_reformat, proteome_df_original, n = 100, protein_header = "PG.ProteinGroups", 
                                          protein_subset_to_use = None, use_only_complete_columns = False):
    proteome_df_aq_reformat = proteome_df_aq_reformat.set_index(["protein", "quant_id"]).replace(0, np.nan)
    if use_only_complete_columns:
        proteome_df_aq_reformat = proteome_df_aq_reformat.dropna()
    
    proteome_df_aq_reformat = proteome_df_aq_reformat.reset_index()

    set_of_proteins = set(proteome_df_aq_reformat["protein"].unique())
    if protein_subset_to_use is not None:
        set_of_proteins = protein_subset_to_use.intersection(set_of_proteins)
    
    return np.random.choice(list(set_of_proteins), n)