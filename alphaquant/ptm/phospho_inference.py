import os
import pathlib
import pandas as pd
import alphaquant.cluster.outlier_scoring as aqoutlier

def get_inferred_phospho_peptides(results_dir, cond1, cond2):
    outlier_handler = aqoutlier.OutlierHandler(results_dir = results_dir, cond1 = cond1, cond2 = cond2)
    clusterdiff_list = outlier_handler.get_diffclust_overview_list()
    predicted_phosphoprone_sequences = load_dl_predicted_phosphoprone_sequences()
    inferred_phospho_peptides = get_regulation_inferred_phosphoprone_peptides(predicted_phosphoprone_sequences, clusterdiff_list)
    return inferred_phospho_peptides

def load_dl_predicted_phosphoprone_sequences(organism = "human"):
    organism_map = {"human": "human_uniprot_reviewed_phos_prob.tsv"}
    database_path = os.path.join(pathlib.Path(__file__).parent.absolute(), "..","resources","phosphopred_databases", organism_map[organism])
    df_phospho_predlib = pd.read_csv(database_path, sep='\t')
    df_phospho_predlib["sequence"] = [f"SEQ_{x}_" for x in df_phospho_predlib["sequence"]]
    return set(df_phospho_predlib[df_phospho_predlib['ptm_prob'] > 0.5]["sequence"])


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
