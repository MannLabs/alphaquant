import os
import pathlib
import pandas as pd
import alphaquant.cluster.outlier_scoring as aqoutlier
import alphaquant.resources.database_loader as aq_resource_dbloader

def get_inferred_phospho_peptides(results_dir, cond1, cond2):
    outlier_handler = aqoutlier.OutlierHandler(results_dir = results_dir, cond1 = cond1, cond2 = cond2)
    clusterdiff_list = outlier_handler.get_diffclust_overview_list()
    predicted_phosphoprone_sequences = aq_resource_dbloader.load_dl_predicted_phosphoprone_sequences()
    inferred_phospho_peptides = get_regulation_inferred_phosphoprone_peptides(predicted_phosphoprone_sequences, clusterdiff_list)
    return inferred_phospho_peptides




def get_regulation_inferred_phosphoprone_peptides(phosphoprone_seqs, clusterdiff_list):
    regulation_inferred_phosphoprone_peptides = []
    for clusterdiff in clusterdiff_list:
        cluster_is_phosphoprone = False
        for seq in clusterdiff.outlier_peptide_names:
            if seq in phosphoprone_seqs:
                cluster_is_phosphoprone = True
                break
        if cluster_is_phosphoprone:
            regulation_inferred_phosphoprone_peptides.extend(clusterdiff.outlier_peptide_names)
    return set(regulation_inferred_phosphoprone_peptides)
