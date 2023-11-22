import pandas as pd

def get_low_correlation_proteoform_df(proteoform_df,max_correlation = 0.7, keep_reference_proteoform = True):
    df_proteoform_lowcorr = proteoform_df[proteoform_df['corr_to_ref'] < max_correlation]

    if keep_reference_proteoform:
        df_proteoform_ref = proteoform_df[proteoform_df["is_reference"]]
        protein_names_lowcorr = df_proteoform_lowcorr['protein'].unique()
        df_proteoform_ref_lowcorr = df_proteoform_ref[df_proteoform_ref['protein'].isin(protein_names_lowcorr)]
        return pd.concat([df_proteoform_lowcorr, df_proteoform_ref_lowcorr])

    return df_proteoform_lowcorr