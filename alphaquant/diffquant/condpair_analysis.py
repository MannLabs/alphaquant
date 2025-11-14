import alphaquant.diffquant.background_distributions as aqbg
import alphaquant.diffquant.diff_analysis as aqdiff
import alphaquant.norm.normalization as aqnorm
import alphaquant.plotting.pairwise as aq_plot_pairwise
import alphaquant.diffquant.diffutils as aqutils
import alphaquant.cluster.cluster_ions as aqclust
import alphaquant.classify.classify_precursors as aq_class_precursors
import alphaquant.cluster.ml_reorder as aq_clust_mlreorder
import alphaquant.tables.diffquant_table as aq_tablewriter_protein
import alphaquant.tables.proteoformtable as aq_tablewriter_proteoform
import alphaquant.tables.misctables as aq_tablewriter_runconfig
import alphaquant.cluster.cluster_utils as aqclust_utils
import alphaquant.cluster.cluster_missingval as aq_clust_missingval
import alphaquant.cluster.outlier_filtering as aq_clust_outlier

import pandas as pd
import numpy as np
import os

from collections import defaultdict

import alphaquant.config.config as aqconfig
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

def analyze_condpair(*,runconfig, condpair):
    """Main workflow orchestration for differential analysis of a condition pair.

    This function coordinates the complete analysis pipeline for comparing two conditions:
    1. Loads and filters data for the two conditions
    2. Performs normalization (within and between conditions)
    3. Creates empirical background distributions
    4. Computes ion-level differential statistics
    5. Builds hierarchical trees and performs clustering to identify proteoforms
    6. Applies machine learning quality scoring (if enabled)
    7. Filters outlier peptides (if enabled)
    8. Generates output tables with FDR correction
    9. Creates visualization plots

    Args:
        runconfig: Configuration object containing all analysis parameters (see run_pipeline docstring)
        condpair: Tuple of (condition1_name, condition2_name) to compare

    Returns:
        tuple: (results_df, peptide_df) - DataFrames with protein and peptide-level results
    """
    LOGGER.info(f"start processeing condpair {condpair}")
    prot2diffions = defaultdict(list) #per default maps any key to empty list
    prot2missingval_diffions = defaultdict(list)
    p2z = {}
    ion2clust = {}
    protnodes = []
    protnodes_missingval = []

    input_df_local = get_unnormed_df_condpair(input_file=runconfig.input_file_reformat, samplemap_df=runconfig.samplemap_df, condpair=condpair, file_has_alphaquant_format = runconfig.file_has_alphaquant_format)
    pep2prot = dict(zip(input_df_local.index, input_df_local['protein']))
    c1_samples, c2_samples = aqutils.get_samples_used_from_samplemap_df(runconfig.samplemap_df, condpair[0], condpair[1])

    try:
        df_c1, df_c2 = get_per_condition_dataframes(c1_samples, c2_samples, input_df_local, min_valid_values=runconfig.min_valid_values, valid_values_filter_mode=runconfig.valid_values_filter_mode, min_valid_values_c1=runconfig.min_valid_values_c1, min_valid_values_c2=runconfig.min_valid_values_c2)
    except Exception as e:
        LOGGER.info(e)
        return

    df_c1_normed, df_c2_normed = aqnorm.normalize_if_specified(df_c1 = df_c1, df_c2 = df_c2, c1_samples = c1_samples, c2_samples = c2_samples, normalize_within_conds = runconfig.normalize, normalize_between_conds = runconfig.normalize,
    runtime_plots = runconfig.runtime_plots, protein_subset_for_normalization_file=runconfig.protein_subset_for_normalization_file, pep2prot = pep2prot)#, "./test_data/normed_intensities.tsv")

    if runconfig.results_dir != None:
        write_out_normed_df(df_c1_normed, df_c2_normed, pep2prot, runconfig.results_dir, condpair)
    normed_c1 = aqbg.ConditionBackgrounds(df_c1_normed, p2z)
    normed_c2 = aqbg.ConditionBackgrounds(df_c2_normed, p2z)

    ions_to_check = normed_c1.ion2nonNanvals.keys() & normed_c2.ion2nonNanvals.keys()
    ions_to_check = sorted(ions_to_check)

    bgpair2diffDist = {}
    deedpair2doublediffdist = {}
    count_ions=0
    all_diffions = []
    for ion in ions_to_check:
        vals1 = normed_c1.ion2nonNanvals.get(ion)
        vals2 = normed_c2.ion2nonNanvals.get(ion)
        bg1 = normed_c1.ion2background.get(ion)
        bg2 = normed_c2.ion2background.get(ion)
        diffDist = aqbg.get_subtracted_bg(bgpair2diffDist, bg1, bg2, p2z)
        if runconfig.ion_test_method == 'ttest':
            diffIon = aqdiff.DifferentialIonTTest(vals1, vals2, ion, p2z, runconfig.outlier_correction)
        else:
            diffIon = aqdiff.DifferentialIon(vals1, vals2, diffDist, ion, runconfig.outlier_correction)
        protein = pep2prot.get(ion)
        if diffIon.usable:
            prot2diffions[protein].append(diffIon)
            all_diffions.append(diffIon)
        else:
            prot2missingval_diffions[protein].append(diffIon)



        if count_ions%2000==0:
            LOGGER.info(f"checked {count_ions} of {len(ions_to_check)} ions")

        count_ions+=1

    count_prots = 0
    for prot in prot2diffions.keys():
        ions = prot2diffions.get(prot)
        if len(ions)<runconfig.min_num_ions:
            continue



        clustered_prot_node = aqclust.get_scored_clusterselected_ions(prot, ions, normed_c1, normed_c2, bgpair2diffDist, p2z, deedpair2doublediffdist,
                                                                        pval_threshold_basis = runconfig.cluster_threshold_pval, fcfc_threshold = runconfig.cluster_threshold_fcfc,
                                                                        take_median_ion=runconfig.take_median_ion, fcdiff_cutoff_clustermerge= runconfig.fcdiff_cutoff_clustermerge,
                                                                        fragment_outlier_filtering=runconfig.fragment_outlier_filtering)
        protnodes.append(clustered_prot_node)

        if count_prots%100==0:
            LOGGER.info(f"checked {count_prots} of {len(prot2diffions.keys())} prots")

        count_prots+=1

    if len(prot2missingval_diffions.keys())>0:
        LOGGER.info(f"start analysis of proteins w. completely missing values")

        for prot in prot2missingval_diffions.keys():

            if prot in prot2diffions.keys(): #only do the missingval analysis if the protein was not analyzed the intensity-based way
                continue
            ions = prot2missingval_diffions.get(prot)
            protnode_missingval = aq_clust_missingval.create_protnode_from_missingval_ions(gene_name=prot,diffions=ions, normed_c1=normed_c1, normed_c2=normed_c2)
            if (protnode_missingval.c1_has_values) and (protnode_missingval.c2_has_values): #one of the conditions has to be missing, otherwise it means that there was e.g. one fragment ion with values in c1 and other fragment ions with values in c2
                continue
            protnodes_missingval.append(protnode_missingval)

        LOGGER.info(f"finished missing value analysis")

    if runconfig.use_ml:
        ml_performance_dict = {}

        #aq_class_stacked_frag.assign_predictability_scores_stacked(protein_nodes= protnodes, acquisition_info_df=None,results_dir=runconfig.results_dir, name = aqutils.get_condpairname(condpair)+"_fragions",
         #                           min_num_fragions=5, replace_nans=True, performance_metrics=ml_performance_dict, plot_predictor_performance=True)
        ml_successfull =aq_class_precursors.assign_predictability_scores(protein_nodes= protnodes, results_dir=runconfig.results_dir, name = aqutils.get_condpairname(condpair), ml_info_file=runconfig.ml_input_file,
                                        samples_used =c1_samples + c2_samples, min_num_precursors=3, prot_fc_cutoff=0, replace_nans=True, performance_metrics=ml_performance_dict, plot_predictor_performance=runconfig.runtime_plots)


        if ml_successfull and (ml_performance_dict["r2_score"] >0.05): #only use the ml score if it is meaningful
            aq_clust_mlreorder.update_nodes_w_ml_score(protnodes)
            LOGGER.info(f"ML based quality score above quality threshold and added to the nodes.")
            runconfig.ml_based_quality_score = True
        else:
            LOGGER.info(f"ML based quality score below quality threshold and not added to the nodes.")
            runconfig.ml_based_quality_score = False

    if runconfig.peptide_outlier_filtering:
        aq_clust_outlier.apply_peptide_outlier_filtering(protnodes)

    protnodes_combined = protnodes + protnodes_missingval
    condpair_node = aqclust_utils.get_condpair_node(protnodes_combined, condpair)
    condpair_node.fraction_missingval = len(prot2missingval_diffions.keys())/(len(prot2diffions.keys())+len(prot2missingval_diffions.keys()))
    res_df, pep_df = write_out_tables(condpair_node, runconfig)

    LOGGER.info(f"condition pair {condpair} finished!")

    return res_df, pep_df

import alphaquant.diffquant.diffutils as aqutils
def get_unnormed_df_condpair(input_file:str, samplemap_df:pd.DataFrame, condpair:str, file_has_alphaquant_format: bool) -> pd.DataFrame:


    samples_c1, samples_c2 = aqutils.get_samples_used_from_samplemap_df(samplemap_df=samplemap_df, cond1 = condpair[0], cond2 = condpair[1])
    used_samples = samples_c1+samples_c2
    unnormed_df = aqutils.import_data(input_file, samples_subset=used_samples, file_has_alphaquant_format = file_has_alphaquant_format)
    unnormed_df, _ = aqutils.prepare_loaded_tables(unnormed_df, samplemap_df)
    return unnormed_df



def write_out_normed_df(normed_df_1, normed_df_2, pep2prot, results_dir, condpair):
    merged_df = normed_df_1.merge(normed_df_2, left_index = True, right_index = True)
    merged_df = 2**merged_df
    merged_df = merged_df.replace(np.nan, 0)
    merged_df["protein"] = list(map(lambda x : pep2prot.get(x),merged_df.index))
    if not os.path.exists(f"{results_dir}/"):
        os.makedirs(f"{results_dir}/")
    merged_df.to_csv(f"{results_dir}/{aqutils.get_condpairname(condpair)}.normed.tsv", sep = "\t")


def get_per_condition_dataframes(samples_c1, samples_c2, unnormed_df, min_valid_values, valid_values_filter_mode, min_valid_values_c1, min_valid_values_c2):

    min_samples = min(len(samples_c1), len(samples_c2))

    if min_samples<2:
        raise Exception(f"condpair has not enough samples: c1:{len(samples_c1)} c2: {len(samples_c2)}, skipping")

    if valid_values_filter_mode == "either":
        min_valid_values = np.min([get_min_valid_values_for_cond(samples_c1, min_valid_values), get_min_valid_values_for_cond(samples_c2, min_valid_values)])
        passes_min_valid_values_c1 = unnormed_df.loc[:, samples_c1].notna().sum(axis=1) >= min_valid_values
        passes_min_valid_values_c2 = unnormed_df.loc[:, samples_c2].notna().sum(axis=1) >= min_valid_values
        passes_min_valid_values = passes_min_valid_values_c1 | passes_min_valid_values_c2
        unnormed_df = unnormed_df[passes_min_valid_values]
        df_c1 = unnormed_df.loc[:, samples_c1]
        df_c2 = unnormed_df.loc[:, samples_c2]

    elif valid_values_filter_mode == "both":
        min_valid_values_c1 = get_min_valid_values_for_cond(samples_c1, min_valid_values)
        min_valid_values_c2 = get_min_valid_values_for_cond(samples_c2, min_valid_values)
        df_c1 = unnormed_df.loc[:, samples_c1].dropna(thresh=min_valid_values_c1, axis=0)
        df_c2 = unnormed_df.loc[:, samples_c2].dropna(thresh=min_valid_values_c2, axis=0)

    elif valid_values_filter_mode == "per_condition":
        min_valid_values_c1 = get_min_valid_values_for_cond(samples_c1, min_valid_values_c1)
        min_valid_values_c2 = get_min_valid_values_for_cond(samples_c2, min_valid_values_c2)
        df_c1 = unnormed_df.loc[:, samples_c1].dropna(thresh=min_valid_values_c1, axis=0)
        df_c2 = unnormed_df.loc[:, samples_c2].dropna(thresh=min_valid_values_c2, axis=0)
    else:
        raise Exception(f"invalid value set for the variable valid_values_filter_mode: {valid_values_filter_mode}, please ensure that is set to: 'either', 'both' or 'per_condition'")

    if (len(df_c1.index)<5) | (len(df_c2.index)<5):
        raise Exception(f"condpair has not enough data for processing c1: {len(df_c1.index)} c2: {len(df_c2.index)}, skipping")

    return df_c1, df_c2

def get_min_valid_values_for_cond(c_samples, min_valid_values):
    if min_valid_values is None: #in the case of None, no nans will be allowed
        return None
    num_samples = len(c_samples)
    if num_samples<min_valid_values:
        return num_samples
    else:
        return min_valid_values




def write_out_tables(condpair_node, runconfig):
    condpair = condpair_node.name

    res_df = aq_tablewriter_protein.TableFromNodeCreator(condpair_node, node_type = "gene", min_num_peptides = runconfig.minpep, annotation_file= getattr(runconfig, "annotation_file", None)).results_df
    has_sequence_nodes = check_if_has_sequence_nodes(condpair_node)
    if has_sequence_nodes:
        pep_df = aq_tablewriter_protein.TableFromNodeCreator(condpair_node, node_type = "seq").results_df
    else:
        pep_df = None
    has_precursor_nodes = check_if_has_precursor_nodes(condpair_node)
    if has_precursor_nodes:
        prec_df = aq_tablewriter_protein.TableFromNodeCreator(condpair_node, node_type = "mod_seq_charge").results_df
    else:
        prec_df = None

    has_base_nodes = check_if_has_base_nodes(condpair_node)
    if has_base_nodes and runconfig.write_base_ions:
        base_df = aq_tablewriter_protein.TableFromNodeCreator(condpair_node, node_type = "base").results_df
    else:
        base_df = None


    if runconfig.runtime_plots:
        aq_plot_pairwise.volcano_plot(res_df, fdr_cutoff= runconfig.volcano_fdr, log2fc_cutoff = runconfig.volcano_fcthresh)
        if has_sequence_nodes:
            aq_plot_pairwise.volcano_plot(pep_df,fdr_cutoff = runconfig.volcano_fdr, log2fc_cutoff = runconfig.volcano_fcthresh)

    if runconfig.results_dir!=None:
        if runconfig.write_out_results_tree:
            aqclust_utils.export_condpairtree_to_json(condpair_node, results_dir = runconfig.results_dir)
        # Write proteoform table defensively; skip if errors occur
        try:
            proteoform_df = aq_tablewriter_proteoform.ProteoFormTableCreator(
                condpair_tree=condpair_node, organism=runconfig.organism
            ).proteoform_df
            proteoform_df.to_csv(
                f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.proteoforms.tsv",
                sep='\t', index=False
            )
        except Exception as e:
            LOGGER.warning(f"Skipping proteoform table write due to error: {e}")

        runconfig_df = aq_tablewriter_runconfig.RunConfigTableCreator(runconfig).runconfig_df

        runconfig_df.to_csv(f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.runconfig.tsv", sep='\t', header=False)
        res_df.to_csv(f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.results.tsv", sep = "\t", index=None)
        if has_sequence_nodes:
            pep_df.to_csv(f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.results.seq.tsv", sep = "\t", index=None)

        if has_precursor_nodes:
            prec_df.to_csv(f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.results.prec.tsv", sep = "\t", index=None)

        if base_df is not None:
            base_df.to_csv(f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.results.base.tsv", sep = "\t", index=None)

    return res_df, pep_df

def check_if_has_sequence_nodes(condpair_node):
    return condpair_node.children[0].children[0].type == "seq"

def check_if_has_precursor_nodes(condpair_node):
    try:
        return condpair_node.children[0].children[0].children[0].children[0].type == "mod_seq_charge"
    except:
        return False

def check_if_has_base_nodes(condpair_node):
    try:
        # Check if we have base nodes (fragments/MS1) at the leaf level
        return condpair_node.children[0].leaves[0].type == "base"
    except:
        return False
