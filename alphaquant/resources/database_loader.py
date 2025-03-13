import os
import pathlib
import pandas as pd
import shutil
import tempfile
import alphabase.tools.data_downloader as ab_downloader
import alphaquant.config.config as aqconfig
import alphaquant.config.variables as aq_variables

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

def _get_path_to_database(database_name, organism):
    """Get the path to a reference database file, which contains mapping to reference protein identifiers, downloading it if necessary.

    Args:
        database_name (str): Name of the database file to locate (e.g., 'uniprot_mapping.tsv', 'swissprot_mapping.tsv')
        organism (str): Organism name (e.g., 'human') - used to locate the correct subfolder

    Raises:
        Exception: If downloading the reference databases from the datashare URL fails, typically due to
                  network issues or invalid datashare access
        Exception: If the requested database file is not found in the expected folder structure after
                  download attempt (/<database_folder>/<organism>/<database_name>)

    Returns:
        str: Full path to the requested database file
    """
    database_folder = os.path.join(pathlib.Path(__file__).parent.absolute(), "reference_databases")
    LOGGER.info(f"Checking for reference databases in {database_folder}")
    if not os.path.exists(database_folder):
        LOGGER.info(f"Downloading reference databases to {database_folder}")
        try:
            datashare_downloader = ab_downloader.DataShareDownloader("https://datashare.biochem.mpg.de/s/ezPzeqStEgDD8gg", output_dir=f"{database_folder}/..") #the url points to the shared AlphaX database, subfolder alphaquant/reference_databases
            datashare_downloader.download()
        except Exception as e:
            error_msg = f"Failed to download reference databases: {str(e)}"
            raise Exception(error_msg) from e
    database_path =  os.path.join(database_folder, organism, database_name)
    if not os.path.exists(database_path):
        raise Exception(f"Reference database {database_name} for organism {organism} not found at {database_path}")
    return database_path


def load_dl_predicted_phosphoprone_sequences(organism = "human"):
    """
    Load the database of peptide sequences, each of which has a "phospho-prone" probability as described
    in the AlphaQuant manuscript. Then filter for those with a phosphorylation probability > 0.5.
    and adapt the sequence to the format required by the AlphaQuant pipeline.

    Args:
        organism (str, optional): The organism to load phospho-prone sequences for. Defaults to "human".

    Raises:
        Exception: If the required database folder structure is not found in the datashare or if
                  downloading the phosphopred databases from the datashare URL fails. This can happen
                  if the datashare is not accessible or if the folder structure doesn't match the expected path.

    Returns:
        set: A set of peptide sequences that have a phosphorylation probability greater than 0.5,
             formatted with "SEQ_" prefix and "_" suffix for AlphaQuant pipeline compatibility.
    """
    organism_map = {"human": "human_uniprot_reviewed_phos_prob.tsv"}
    database_folder = os.path.join(pathlib.Path(__file__).parent.absolute(), "..","resources","phosphopred_databases")

    LOGGER.info(f"Checking for phosphopred databases in {database_folder}")
    if not os.path.exists(database_folder):
        LOGGER.info(f"Downloading phosphopred databases to {database_folder}")
        try:
            datashare_downloader = ab_downloader.DataShareDownloader("https://datashare.biochem.mpg.de/s/stH9pmNe6O9CRHG", output_dir=f"{database_folder}/..") #the url points to the shared AlphaX database, subfolder alphaquant/phosphopred_databases
            datashare_downloader.download()
        except Exception as e:
            error_msg = f"Failed to download phosphopred databases: {str(e)}"
            raise Exception(error_msg) from e
    database_path = os.path.join(database_folder, organism_map[organism])

    df_phospho_predlib = pd.read_csv(database_path, sep='\t')
    df_phospho_predlib["sequence"] = [f"{aq_variables.SEQ}_{x}_" for x in df_phospho_predlib["sequence"]]
    return set(df_phospho_predlib[df_phospho_predlib['ptm_prob'] > 0.5]["sequence"])
