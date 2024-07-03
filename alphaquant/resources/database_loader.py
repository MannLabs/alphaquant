import os
import pathlib
import pandas as pd

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
    swissprot_df = pd.read_csv(swissprot_file, sep = '\t', usecols=["Gene Names", 'Entry'])
    gene_names = swissprot_df["Gene Names"].astype(str).tolist()
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
    database_path =  os.path.join(pathlib.Path(__file__).parent.absolute(), "reference_databases", organism, database_name)
    return database_path




