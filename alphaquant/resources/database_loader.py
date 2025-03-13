import os
import pathlib
import pandas as pd
import shutil
import tempfile
import alphabase.tools.data_downloader as ab_downloader
import alphaquant.config.config as aqconfig
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)


def get_genename2sequence_dict( organism = "human"):
    swissprot_file = get_swissprot_path(organism)

    swissprot_df = pd.read_csv(swissprot_file, sep = '\t', usecols=["Gene Names", 'Sequence'])
    gene_names = swissprot_df["Gene Names"].astype(str).tolist()

    sequences = swissprot_df['Sequence'].astype(str).tolist()

    gene2sequence_dict = {}

    for gene_group, sequence in zip(gene_names, sequences):
        for gene in gene_group.split(" "):
            gene2sequence_dict[gene] = sequence

    return gene2sequence_dict

def get_swissprot2sequence_dict( organism = "human"):
    swissprot_file = get_swissprot_path(organism)
    swissprot_df = pd.read_csv(swissprot_file, sep = '\t', usecols=['Entry', 'Sequence'])
    swissprot_ids = swissprot_df['Entry'].astype(str).tolist()
    sequences = swissprot_df['Sequence'].astype(str).tolist()

    swissprot2sequence_dict = dict(zip(swissprot_ids, sequences))
    return swissprot2sequence_dict

def get_uniprot2sequence_dict( organism = "human"):
    swissprot_file = get_swissprot_path(organism)
    swissprot_df = pd.read_csv(swissprot_file, sep = '\t', usecols=['Entry', 'Sequence'])
    swissprot_ids = swissprot_df['Entry'].astype(str).tolist()
    sequences = swissprot_df['Sequence'].astype(str).tolist()

    swissprot2sequence_dict = dict(zip(swissprot_ids, sequences))
    return swissprot2sequence_dict

def get_genename2swissprot_dict( organism = "human"):
    swissprot_file = get_swissprot_path(organism)
    swissprot_df = pd.read_csv(swissprot_file, sep = '\t', usecols=['Gene Names', 'Entry'])
    gene_names = swissprot_df['Gene Names'].astype(str).tolist()

    swissprot_ids = swissprot_df['Entry'].astype(str).tolist()

    gene2swissprot_dict = {}

    for gene_group, entry in zip(gene_names, swissprot_ids):
        for gene in gene_group.split(" "):
            gene2swissprot_dict[gene] = entry
    return gene2swissprot_dict



def get_uniprot_path( organism= "human"):
    return _get_path_to_database("uniprot_mapping.tsv",organism)

def get_swissprot_path( organism = "human"):
    return _get_path_to_database("swissprot_mapping.tsv",organism)

def _get_path_to_database( database_name, organism):
    database_folder = os.path.join(pathlib.Path(__file__).parent.absolute(), "reference_databases")
    LOGGER.info(f"Checking for reference databases in {database_folder}")
    if not os.path.exists(database_folder):
        LOGGER.info(f"Downloading reference databases to {database_folder}")
        try:
            dsd = ab_downloader.DataShareDownloader("https://datashare.biochem.mpg.de/s/ezPzeqStEgDD8gg", output_dir=f"{database_folder}/..") #the url points to the shared AlphaX database, subfolder alphaquant/reference_databases
            dsd.download()
        except Exception as e:
            LOGGER.error(f"Failed to download reference databases: {str(e)}")
            raise Exception(f"Failed to download reference databases: {str(e)}") from e
    database_path =  os.path.join(database_folder, organism, database_name)
    if not os.path.exists(database_path):
        raise Exception(f"Reference database {database_name} for organism {organism} not found at {database_path}")
    return database_path


def load_dl_predicted_phosphoprone_sequences(organism = "human"):
    organism_map = {"human": "human_uniprot_reviewed_phos_prob.tsv"}
    database_folder = os.path.join(pathlib.Path(__file__).parent.absolute(), "..","resources","phosphopred_databases")

    LOGGER.info(f"Checking for phosphopred databases in {database_folder}")
    if not os.path.exists(database_folder):
        LOGGER.info(f"Downloading phosphopred databases to {database_folder}")
        try:
            dsd = ab_downloader.DataShareDownloader("https://datashare.biochem.mpg.de/s/stH9pmNe6O9CRHG", output_dir=f"{database_folder}/..") #the url points to the shared AlphaX database, subfolder alphaquant/phosphopred_databases
            dsd.download()
        except Exception as e:
            LOGGER.error(f"Failed to download phosphopred databases: {str(e)}")
            raise Exception(f"Failed to download phosphopred databases: {str(e)}") from e
    database_path = os.path.join(database_folder, organism_map[organism])

    df_phospho_predlib = pd.read_csv(database_path, sep='\t')
    df_phospho_predlib["sequence"] = [f"SEQ_{x}_" for x in df_phospho_predlib["sequence"]]
    return set(df_phospho_predlib[df_phospho_predlib['ptm_prob'] > 0.5]["sequence"])
