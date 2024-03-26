import os

os.chdir("analyses_data/quicktests/mixed_species_medium_size")

INPUT_FILE = "20210210_154121_S209-S-1-240min_Report_quicktest_shortened_mediumsize.tsv"
HOUSEKEEPING_PROTEINS = "housekeeping_proteins.tsv"
SAMPLEMAP = "samplemap.tsv"
RESULTS_DIR = "results"
SHARED_PEPTIDES_BETWEEN_SPECIES_FILE = "../../databases/intersecting_peptides_human_yeast_cael_ecoli.tsv"

import alphaquant.run_pipeline as run_pipeline

run_pipeline.run_pipeline(input_file=INPUT_FILE, samplemap_file=SAMPLEMAP, results_dir=RESULTS_DIR, runtime_plots=True, minrep_either= 2, protein_subset_for_normalization_file=HOUSEKEEPING_PROTEINS, take_median_ion= True,
                           annotation_columns=["PG.Genes", "PG.Organisms"], input_type_to_use= "spectronaut_fragion_isotopes_protein", peptides_to_exclude_file=SHARED_PEPTIDES_BETWEEN_SPECIES_FILE)