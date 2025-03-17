import os
import re
import pathlib
import pandas as pd
from itertools import combinations
from typing import Optional, List, Tuple, Union
import pandas as pd

import alphaquant.diffquant.diffutils as aq_diffquant_utils
import alphaquant.utils.utils as aq_utils

import alphaquant.ptm.ptmsite_mapping as aqptm
import multiprocess
import alphaquant.config.variables as aqvariables
import alphabase.quantification.quant_reader.config_dict_loader as config_dict_loader
config_dict_loader.INTABLE_CONFIG = os.path.join(pathlib.Path(__file__).parent.absolute(), "./config/quant_reader_config.yaml")

import alphaquant.classify.ml_info_table as aq_ml_info_table
import alphabase.quantification.quant_reader.quant_reader_manager as abquantreader
import alphaquant.tables.alphadia_reader as aq_table_alphadiareader
import alphaquant.diffquant.condpair_analysis as aqcondpair
import alphaquant.multicond.median_condition_creation as aqmediancreation
import alphaquant.multicond.median_condition_analysis as aqmediancond
import alphaquant.tables.misctables as aq_tablewriter_misc

import logging
import shutil
LOGGER = logging.getLogger(__name__)



class ConfigOfRunPipeline:
    """Stores all the parameters given to the "run_pipeline" function"""
    def __init__(self, locals):
        for k, v in locals.items():
            setattr(self, k, v)


def run_pipeline(input_file: str,
                samplemap_file: Optional[str] = None,
                samplemap_df: Optional[pd.DataFrame] = None,
                modification_type: Optional[str] = None,
                input_type_to_use: Optional[str] = None,
                results_dir: str = "./results",
                multicond_median_analysis: bool = False,
                condpairs_list: Optional[List[Tuple[str, str]]] = None,
                file_has_alphaquant_format: bool = False,
                minrep_both: int = 2,
                minrep_either: Optional[int] = None,
                minrep_c1: Optional[int] = None,
                minrep_c2: Optional[int] = None,
                min_num_ions: int = 1,
                minpep: int = 1,
                organism: Optional[str] = None,
                cluster_threshold_pval: float = 0.001,
                cluster_threshold_fcfc: float = 0,
                fcdiff_cutoff_clustermerge = 0.5,
                use_ml: bool = True,
                take_median_ion: bool = True,
                perform_ptm_mapping: bool = False,
                perform_phospho_inference: bool = False,
                outlier_correction: bool = True,
                normalize: bool = True,
                use_iontree_if_possible: bool = True,
                write_out_results_tree: bool = True,
                use_multiprocessing: bool = False,
                runtime_plots: bool = True,
                volcano_fdr: float = 0.05,
                volcano_fcthresh: float = 0.5,
                annotation_columns: Optional[List[str]] = None,
                protein_subset_for_normalization_file: Optional[str] = None,
                protnorm_peptides: bool = True,
                peptides_to_exclude_file: Optional[str] = None,
                reset_progress_folder: bool = False) -> None:
    """Run differential analyses following the AlphaQuant pipeline. This function processes proteomics data through multiple steps including
        preprocessing, if applicable PTM site mapping, if applicable median condition creation, normalization, statistical testing, visualizations
        and writing of results tables.

    Args:
    input_file (str): Proteomics results file from a search engine (e.g., DIA-NN, Spectronaut, AlphaDIA). Required.
    samplemap_file (str): Path to TSV file mapping samples to conditions. Required if samplemap_df not provided.
    samplemap_df (pandas.DataFrame): DataFrame mapping samples to conditions. Required if samplemap_file not provided.
    modification_type (str): For PTM analysis, specifies the modification type (e.g., '[Phospho (STY)]' for Spectronaut). Required if perform_ptm_mapping is True.
    input_type_to_use (str): Type of quantitative information to extract from results table. See intable.config for options.
    results_dir (str): Output directory for analysis results. Defaults to "./results".
    multicond_median_analysis (bool): Whether to compare all conditions to a median condition. Defaults to False.
    condpairs_list (list): Specific condition pairs to compare. If None, performs all pairwise comparisons.
    file_has_alphaquant_format (bool): Whether the input file is already in AlphaQuant matrix format. Defaults to False.
    minrep_both (int): Minimum replicate count required in both conditions. Defaults to 2.
    minrep_either (int): Minimum replicate count required in either condition.
    minrep_c1 (int): Minimum replicate count required in condition 1.
    minrep_c2 (int): Minimum replicate count required in condition 2.
    min_num_ions (int): Minimum number of ions required per peptide. Defaults to 1.
    minpep (int): Minimum number of peptides required per protein. Defaults to 1.
    organism (str): Organism name for PTM mapping (e.g., 'human', 'mouse'). Required if perform_ptm_mapping is True.
    cluster_threshold_pval (float): P-value threshold for statistical clustering. Defaults to 0.001.
    cluster_threshold_fcfc (float): Fold change threshold for clustering. Defaults to 0.
    fcdiff_cutoff_clustermerge (float): Fold change difference cutoff for merging peptide clusters. Defaults to 0.5.
    use_ml (bool): Enable machine learning analysis. Defaults to True.
    take_median_ion (bool): Use median-centered fragment ions for peptide comparisons. Defaults to True.
    perform_ptm_mapping (bool): Enable PTM site mapping analysis. Defaults to False.
    perform_phospho_inference (bool): Enable phosphorylation-prone region annotation. Defaults to False.
    outlier_correction (bool): Enable outlier correction in differential testing. Defaults to True.
    normalize (bool): Enable sample and condition normalization. Defaults to True.
    use_iontree_if_possible (bool): Use ion tree structure when available. Defaults to True.
    write_out_results_tree (bool): Write results in hierarchical tree format. Defaults to True.
    use_multiprocessing (bool): Enable parallel processing. Defaults to False.
    runtime_plots (bool): Generate diagnostic plots including volcanos. Defaults to True.
    volcano_fdr (float): FDR cutoff for volcano plot significance. Defaults to 0.05.
    volcano_fcthresh (float): Fold change threshold for volcano plot significance. Defaults to 0.5.
    annotation_columns (list): Additional columns to include in output tables.
    protein_subset_for_normalization_file (str): File specifying proteins to use for normalization.
    protnorm_peptides (bool): Enable protein-level peptide normalization. Defaults to True.
    peptides_to_exclude_file (str): File listing peptides to exclude (e.g., shared between species).
    reset_progress_folder (bool): Clear and recreate the progress folder. Defaults to False.
    """
    LOGGER.info("Starting AlphaQuant")
    input_file_original = input_file
    check_input_consistency(input_file_original, samplemap_file, samplemap_df)
    create_progress_folder_if_applicable(input_file_original, reset_progress_folder)

    if samplemap_df is None:
        samplemap_df = aq_diffquant_utils.load_samplemap(samplemap_file)

    input_type, config_dict, _ = config_dict_loader.get_input_type_and_config_dict(input_file_original, input_type_to_use)
    annotation_file = load_annotation_file(input_file_original, input_type, annotation_columns)
    use_ml = check_if_table_supports_ml(config_dict)

    if perform_ptm_mapping:
        if modification_type is None:
            raise Exception("modification_type is None, but perform_ptm_mapping is True. Please set perform_ptm_mapping to False or specify modification_type.")
        input_file_reformat = load_ptm_input_file(input_file = input_file_original, input_type_to_use = "spectronaut_ptm_fragion", results_dir = results_dir, samplemap_df = samplemap_df, modification_type = modification_type, organism = organism)
        if use_ml:
            ml_input_file = load_ml_info_file(input_file_original, input_type, modification_type)

    elif "fragment_precursorfiltered.matrix" in input_file_original:
        alphadia_tableprocessor = aq_table_alphadiareader.AlphaDIAFragTableProcessor(input_file_original)
        input_file_reformat = alphadia_tableprocessor.input_file_reformat
        if use_ml:
            ml_input_file = alphadia_tableprocessor.ml_info_file
    else:
        input_file_reformat = load_input_file(input_file_original, input_type)
        if use_ml:
            ml_input_file = load_ml_info_file(input_file_original, input_type)

    if peptides_to_exclude_file is not None:
        remove_peptides_to_exclude_from_input_file(input_file_reformat, peptides_to_exclude_file)

    if multicond_median_analysis:
        condpairs_list = aqmediancreation.get_all_conds_relative_to_median(samplemap_df)
        median_manager = aqmediancreation.MedianConditionManager(input_file_reformat, samplemap_file) #writes median condition to input file and samplemap file and overwrites the formatted input and samplemap file
        input_file_reformat = median_manager.input_filename_adapted
        samplemap_df = median_manager.samplemap_df_extended
        del median_manager #delete the object as it needs not be in the runconfig

    aqvariables.determine_variables(input_file_reformat, input_type)

    #use runconfig object to store the parameters
    runconfig = ConfigOfRunPipeline(locals()) #all the parameters given into the function are transfered to the runconfig object!

    #store method parameters for reproducibility
    aq_diffquant_utils.remove_old_method_parameters_file_if_exists(results_dir)
    aq_diffquant_utils.store_method_parameters(locals(), results_dir)


    if condpairs_list == None:
        conds = samplemap_df["condition"].unique()
        conds = sorted(conds)
        condpairs_list = combinations(conds, 2)

    num_cores = get_num_cores_to_use(use_multiprocessing)

    if num_cores == 1:
        run_analysis_singleprocess(condpair_combinations=condpairs_list, runconfig=runconfig)

    else:
        run_analysis_multiprocess(condpair_combinations=condpairs_list, runconfig=runconfig, num_cores=num_cores)

    if multicond_median_analysis:
        aqmediancond.analyze_and_write_median_condition_results(results_dir)


def check_input_consistency(input_file, samplemap_file, samplemap_df):
    if input_file is None:
        raise Exception("no input file!")
    if samplemap_file is None and samplemap_df is None:
        raise Exception("Samplemap is missing!")
    return True

def create_progress_folder_if_applicable(input_file, reset_progress_folder):
    progress_folder = os.path.join(os.path.dirname(input_file), "progress")
    if not os.path.exists(progress_folder):
        os.makedirs(progress_folder)
    elif reset_progress_folder:
        shutil.rmtree(progress_folder)
        os.makedirs(progress_folder)


def load_ptm_input_file(input_file, input_type_to_use, results_dir, samplemap_df, modification_type, organism):
    reformatted_input_filename = aq_utils.get_progress_folder_filename(input_file, f".ptmsite_mapped.tsv.{input_type_to_use}.aq_reformat.tsv", remove_extension=True)
    if os.path.exists(reformatted_input_filename):#in case there already is a reformatted file, we don't need to reformat it again
        LOGGER.info(f"Reformatted input file already exists. Using reformatted file of type {input_type_to_use}")
        return reformatted_input_filename
    else:
        ptm_mapped_file = write_ptm_mapped_input(input_file, results_dir, samplemap_df, modification_type, organism)
        return load_input_file(ptm_mapped_file, input_type_to_use)

def write_ptm_mapped_input(input_file, results_dir, samplemap_df, modification_type, organism = "human"):
    try:
        aqptm.assign_dataset_inmemory(input_file = input_file, results_dir=results_dir, samplemap_df=samplemap_df, modification_type=modification_type, organism=organism)
    except Exception as e:
        LOGGER.error(f"PTM mapping in memory failed with error: {e}. Trying out-of-core approach with dask.")
        aqptm.assign_dataset_chunkwise(input_file = input_file, results_dir=results_dir, samplemap_df=samplemap_df, modification_type=modification_type, organism=organism)
    mapped_df = pd.read_csv(f"{results_dir}/ptm_ids.tsv", sep = "\t")
    ptm_mapped_file = aqptm.merge_ptmsite_mappings_write_table(input_file, mapped_df, modification_type)
    return ptm_mapped_file


def load_input_file(input_file, input_type):
    reformatted_input_filename = aq_utils.get_progress_folder_filename(input_file, f".{input_type}.aq_reformat.tsv", remove_extension=False)
    if os.path.exists(reformatted_input_filename):#in case there already is a reformatted file, we don't need to reformat it again
        LOGGER.info(f"Reformatted input file already exists. Using reformatted file of type {input_type}")
        return reformatted_input_filename
    else:
        reformatted_input_file_initial = abquantreader.reformat_and_save_input_file(input_file, input_type_to_use = input_type, use_alphaquant_format=True)
        shutil.move(reformatted_input_file_initial, reformatted_input_filename)

    return reformatted_input_filename



def load_annotation_file(input_file, input_type, annotation_columns):
    annotation_filename = aq_utils.get_progress_folder_filename(input_file, f".annotation.tsv")
    if os.path.exists(annotation_filename):#in case there already is a reformatted file, we don't need to reformat it again
        LOGGER.info(f"Annotation file already exists. Using annotation file of type {input_type}")
        return annotation_filename
    else:
        return aq_tablewriter_misc.AnnotationFileCreator(input_file, input_type, annotation_columns).annotation_filename

def check_if_table_supports_ml(config_dict):
    return config_dict["format"] == "longtable"

def load_ml_info_file(input_file, input_type, modification_type = None):
    ml_info_filename = aq_utils.get_progress_folder_filename(input_file, f".ml_info_table.tsv")
    if os.path.exists(ml_info_filename):#in case there already is a reformatted file, we don't need to reformat it again
        LOGGER.info(f"ML info file already exists. Using ML info file of type {input_type}")
        return ml_info_filename
    else:
        return aq_ml_info_table.MLInfoTableCreator(input_file, input_type, modification_type).ml_info_filename


def remove_peptides_to_exclude_from_input_file(input_file, peptides_to_exclude_file):
    df_input = pd.read_csv(input_file, sep = "\t")
    peptides_to_exclude = set(pd.read_csv(peptides_to_exclude_file, sep = "\t")["peptide"].tolist())
    pattern = r"SEQ_([A-Za-z0-9]+)_"
    try:
        df_input["peptide"] = [re.search(pattern, peptide).group(1) for peptide in df_input[aqvariables.QUANT_ID]]
    except:
        raise Exception("parsing of peptide sequence from QUANT_ID failed. The QUANT_ID column should contain the peptide sequence in the format SEQ_<peptide>_")

    not_in_peptides_to_exclude = ~df_input["peptide"].isin(peptides_to_exclude)
    df_input = df_input[not_in_peptides_to_exclude]
    df_input = df_input.drop(columns = ["peptide"])
    df_input.to_csv(input_file, sep = "\t", index = False)
    num_removed = len(not_in_peptides_to_exclude) - len(df_input.index)
    LOGGER.info(f"Excluded {num_removed} shared-species entries from input file")



def get_num_cores_to_use(use_multiprocessing):
    num_cores = multiprocess.cpu_count() if use_multiprocessing else 1
    return min(num_cores, 10)

def run_analysis_singleprocess(condpair_combinations, runconfig):

    for condpair in condpair_combinations:
        aqcondpair.analyze_condpair(runconfig=runconfig, condpair=condpair)

def run_analysis_multiprocess(condpair_combinations, runconfig, num_cores):

    with multiprocess.Pool(num_cores) as pool:

        pool.map(lambda condpair :

        aqcondpair.analyze_condpair(runconfig= runconfig, condpair = condpair)

        ,condpair_combinations)

