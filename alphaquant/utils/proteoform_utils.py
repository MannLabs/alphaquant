import pandas as pd

def filter_proteoform_df(proteoform_df, min_num_peptides=1, quality_score_cutoff=1, fcdiff_cutoff=0.5, likely_phospho = None,keep_reference_proteoform=True):

    proteoform_df_outliers = proteoform_df[(proteoform_df['abs_fcdiff'] > fcdiff_cutoff) & (proteoform_df['quality_score'] <= quality_score_cutoff) & (proteoform_df['num_peptides'] >= min_num_peptides)]
    if likely_phospho is not None:
        proteoform_df_outliers = proteoform_df_outliers[proteoform_df_outliers['likely_phospho'] == likely_phospho]

    if keep_reference_proteoform:
        df_proteoform_ref = proteoform_df[proteoform_df["is_reference"]]
        protein_names_lowcorr = proteoform_df_outliers['protein'].unique()
        df_proteoform_ref_lowcorr = df_proteoform_ref[df_proteoform_ref['protein'].isin(protein_names_lowcorr)]
        return pd.concat([proteoform_df_outliers, df_proteoform_ref_lowcorr])

    return proteoform_df_outliers
