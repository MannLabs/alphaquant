import alphaquant.diffquant.background_distributions as aqbg
import alphaquant.diffquant.diff_analysis as aqdiff
import alphaquant.norm.normalization as aqnorm
import alphaquant.plotting.pairwise as aq_plot_pairwise
import alphaquant.diffquant.diffutils as aqutils
import alphaquant.cluster.cluster_ions as aqclust
import alphaquant.classify.classify_ions as aqclass
import alphaquant.config.variables as aqvariables
import alphaquant.diffquant.tablewriter as aq_diffquant_tablewriter
import anytree

import alphaquant.cluster.cluster_utils as aqclust_utils
import statsmodels.stats.multitest as mt
from time import time
import pandas as pd
import numpy as np
import os

def analyze_condpair(*,runconfig, condpair):
    t_zero = time()
    print(f"start processeing condpair {condpair}")
    prot2diffions = {}
    p2z = {}
    ion2clust = {}
    protnodes = []
    quantified_peptides = []
    quantified_proteins = []


    input_df_local = get_unnormed_df_condpair(input_file=runconfig.input_file, samplemap_df=runconfig.samplemap_df, condpair=condpair, file_has_alphaquant_format = runconfig.file_has_alphaquant_format)
    pep2prot = dict(zip(input_df_local.index, input_df_local['protein']))
    c1_samples, c2_samples = aqutils.get_samples_used_from_samplemap_df(runconfig.samplemap_df, condpair[0], condpair[1])

    try:
        df_c1, df_c2 = get_per_condition_dataframes(c1_samples, c2_samples, input_df_local, runconfig.minrep)
    except Exception as e:
        print(e)
        return

    df_c1_normed, df_c2_normed = aqnorm.normalize_if_specified(df_c1 = df_c1, df_c2 = df_c2, c1_samples = c1_samples, c2_samples = c2_samples, minrep = runconfig.minrep, normalize_within_conds = runconfig.normalize, normalize_between_conds = runconfig.normalize,
    runtime_plots = runconfig.runtime_plots, protein_subset_for_normalization_file=runconfig.protein_subset_for_normalization_file, pep2prot = pep2prot,prenormed_file = runconfig.pre_normed_intensity_file)#, "./test_data/normed_intensities.tsv")

    if runconfig.results_dir != None:
        write_out_normed_df(df_c1_normed,df_c2_normed, pep2prot, runconfig.results_dir, condpair)
    t_normalized = time()
    normed_c1 = aqbg.ConditionBackgrounds(df_c1_normed, p2z)
    normed_c2 = aqbg.ConditionBackgrounds(df_c2_normed, p2z)
    
    t_bgdist_fin = time()
    ions_to_check = normed_c1.ion2nonNanvals.keys() & normed_c2.ion2nonNanvals.keys()
    use_ion_tree = runconfig.use_iontree_if_possible
    bgpair2diffDist = {}
    deedpair2doublediffdist = {}
    count_ions=0
    for ion in ions_to_check:
        t_ion = time()
        vals1 = normed_c1.ion2nonNanvals.get(ion)
        vals2 = normed_c2.ion2nonNanvals.get(ion)
        bg1 = normed_c1.ion2background.get(ion)
        bg2 = normed_c2.ion2background.get(ion)
        diffDist = aqbg.get_subtracted_bg(bgpair2diffDist, bg1, bg2, p2z)
        t_subtract_end = time()
        diffIon = aqdiff.DifferentialIon(vals1, vals2, diffDist, ion, runconfig.outlier_correction)
        t_diffion = time()
        protein = pep2prot.get(ion)
        prot_ions = prot2diffions.get(protein, list())
        prot_ions.append(diffIon)
        prot2diffions[protein] = prot_ions
        quantified_peptide = QuantifiedResult(kwargs = {aqvariables.QUANT_ID : ion, 'pval' : diffIon.p_val, 'log2fc' : diffIon.fc, 'protein' : protein, 'condpair' : aqutils.get_condpairname(condpair)})
        quantified_peptides.append(quantified_peptide)


        if count_ions%2000==0:
            print(f"checked {count_ions} of {len(ions_to_check)} ions")

        count_ions+=1


    count_prots = 0
    for prot in prot2diffions.keys():
        ions = prot2diffions.get(prot)
        if len(ions)<runconfig.min_num_ions:
            continue
        diffprot = aqdiff.DifferentialProtein(prot, ions, runconfig.median_offset, runconfig.dia_fragment_selection)

        if use_ion_tree:
            clustered_root_node = aqclust.get_scored_clusterselected_ions(prot, ions, normed_c1, normed_c2, bgpair2diffDist, p2z, deedpair2doublediffdist, 
                                                                          pval_threshold_basis = runconfig.cluster_threshold_pval, fcfc_threshold = runconfig.cluster_threshold_fcfc, 
                                                                          take_median_ion=runconfig.take_median_ion, fcdiff_cutoff_clustermerge= runconfig.fcdiff_cutoff_clustermerge)
            #print(anytree.RenderTree(clustered_root_node))
            #if not clustered_root_node.is_included:
             #   continue
            protnodes.append(clustered_root_node)
            pval, fc, consistency_score, ions_included = aqclust_utils.get_diffresults_from_clust_root_node(clustered_root_node)
            num_peptides = len(anytree.findall(clustered_root_node, filter_ = lambda x : x.type == 'seq'))
            if num_peptides < runconfig.minpep:
                continue

        else:
            pval, fc, consistency_score, ions_included = diffprot.pval, diffprot.fc, np.nan,diffprot.ions

        if runconfig.get_ion2clust:
            ionclust_protein = aqclust.find_fold_change_clusters_base_ions([[x] for x in ions],normed_c1,normed_c2, bgpair2diffDist,p2z, deedpair2doublediffdist, fc_threshold=runconfig.cluster_threshold_fcfc, pval_threshold_basis=runconfig.cluster_threshold_pval)
            ion2clust.update({x.name:y for x,y in ionclust_protein.items()})

        if count_prots%100==0:
            print(f"checked {count_prots} of {len(prot2diffions.keys())} prots")
        count_prots+=1

        pseudoint1_cond, pseudoint2_cond = aqdiff.calc_pseudo_intensities(ions, normed_c2, diffprot.fc)
        quantified_protein = QuantifiedResult(kwargs={'condpair': aqutils.get_condpairname(condpair), 'protein' : prot, 'fdr' : 1.0, 'pval':pval,'log2fc': fc, 'consistency_score' : consistency_score, 'num_ions' : len(ions_included), 'pseudoint1' : pseudoint1_cond, 'pseudoint2' : pseudoint2_cond})
        quantified_proteins.append(quantified_protein)

    if use_ion_tree:
        if runconfig.use_ml & len(protnodes)>100:
            ml_performance_dict = {}
            ml_successfull = True
            aqclass.assign_predictability_scores(protnodes, runconfig.results_dir, name = aqutils.get_condpairname(condpair), 
                                                 samples_used = c1_samples+ c2_samples,precursor_cutoff=3,
            fc_cutoff=0.5, number_splits=5, plot_predictor_performance=runconfig.runtime_plots, 
            replace_nans=True, performance_metrics=ml_performance_dict, protnorm_peptides=runconfig.protnorm_peptides)


            if (ml_performance_dict["r2_score"] >0.05) and ml_successfull: #only use the ml score if it is meaningful
                aqclust.update_nodes_w_ml_score(protnodes)
                update_quantified_proteins_w_tree_results(quantified_proteins, protnodes)


    add_fdr(quantified_proteins)
    add_fdr(quantified_peptides)

    res_df = get_results_df(quantified_proteins)
    pep_df = get_results_df(quantified_peptides)

    if runconfig.runtime_plots:
        aq_plot_pairwise.volcano_plot(res_df, significance_cutoff = runconfig.volcano_fdr, log2fc_cutoff = runconfig.volcano_fcthresh)
        aq_plot_pairwise.volcano_plot(pep_df,significance_cutoff = runconfig.volcano_fdr, log2fc_cutoff = runconfig.volcano_fcthresh)

    if runconfig.results_dir!=None:

        if runconfig.write_out_results_tree:
            condpair_node = aqclust_utils.get_condpair_node(protnodes, condpair)
            aqclust_utils.export_condpairtree_to_json(condpair_node, results_dir = runconfig.results_dir)
            proteoform_df = aq_diffquant_tablewriter.ProteoFormTableCreator(condpair_tree= condpair_node, organism=runconfig.organism_for_phospho_inference).proteoform_df

            proteoform_df.to_csv(f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.proteoforms.tsv", sep='\t', index=False)

        if runconfig.annotation_file != None: #additional annotations can be added before saving
            annot_df = pd.read_csv(runconfig.annotation_file, sep = "\t")
            intersect_columns = annot_df.columns.intersection(pep_df.columns)
            if(len(intersect_columns)>0):
                print(list(intersect_columns))
                res_df = res_df.merge(annot_df, on=list(intersect_columns), how= 'left')
                pep_df = pep_df.merge(annot_df, on= list(intersect_columns), how = 'left')



        if runconfig.get_ion2clust:
            ion2clust_df = pd.DataFrame(ion2clust.items(), columns=['quant_id', 'cluster'])
            ion2clust_df.to_csv(f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.ion2clust.tsv", sep = "\t", index=None)

        res_df.to_csv(f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.results.tsv", sep = "\t", index=None)
        pep_df.to_csv(f"{runconfig.results_dir}/{aqutils.get_condpairname(condpair)}.results.ions.tsv", sep = "\t", index=None)

    print(f"\ncondition pair {condpair} finished!\n")


    return res_df, pep_df

import alphaquant.diffquant.diffutils as aqutils
def get_unnormed_df_condpair(input_file:str, samplemap_df:pd.DataFrame, condpair:str, file_has_alphaquant_format: bool) -> pd.DataFrame:


    samples_c1, samples_c2 = aqutils.get_samples_used_from_samplemap_df(samplemap_df=samplemap_df, cond1 = condpair[0], cond2 = condpair[1])
    used_samples = samples_c1+samples_c2
    unnormed_df = aqutils.import_data(input_file, samples_subset=used_samples, file_has_alphaquant_format = file_has_alphaquant_format)
    unnormed_df, _ = aqutils.prepare_loaded_tables(unnormed_df, samplemap_df)
    return unnormed_df


# Cell
#helper class to store diffresults
class QuantifiedResult:
    def __init__(self, **kwargs):
        self.propdict = None
        if kwargs:
            self.propdict = kwargs['kwargs']
    def add_property(self, key, value):
        self.propdict[key]  = value
    def add_properties(self, dict):
        self.propdict.update(dict)


def update_quantified_proteins_w_tree_results(quantified_proteins, protnodes):
    prot2fc = {x.name : x.fc for x in protnodes}
    prot2pval = {x.name : x.p_val  for x in protnodes}
    prot2predscore = {x.name : x.predscore for x in protnodes}
    for quantified_protein in quantified_proteins:
        protname = quantified_protein.propdict['protein']
        quantified_protein.propdict['log2fc'] = prot2fc.get(protname)
        quantified_protein.propdict['pval'] = prot2pval.get(protname)
        quantified_protein.propdict['predscore'] = prot2predscore.get(protname)


def add_fdr(quantified_results):
    pvals = [x.propdict["pval"] for x in quantified_results]
    fdrs = mt.multipletests(pvals, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
    for idx in range(len(quantified_results)):
        quantified_results[idx].propdict['fdr'] = fdrs[idx]

def get_results_df(quantfied_results):
    quantified_results_dicts = [x.propdict for x in quantfied_results]
    res_df = pd.DataFrame(quantified_results_dicts)
    return res_df

def write_out_normed_df(normed_df_1, normed_df_2, pep2prot, results_dir, condpair):
    merged_df = normed_df_1.merge(normed_df_2, left_index = True, right_index = True)
    merged_df = 2**merged_df
    merged_df = merged_df.replace(np.nan, 0)
    merged_df["protein"] = list(map(lambda x : pep2prot.get(x),merged_df.index))
    if not os.path.exists(f"{results_dir}/"):
        os.makedirs(f"{results_dir}/")
    merged_df.to_csv(f"{results_dir}/{aqutils.get_condpairname(condpair)}.normed.tsv", sep = "\t")


def get_per_condition_dataframes(samples_c1, samples_c2, unnormed_df, minrep):

    min_samples = min(len(samples_c1), len(samples_c2))

    if min_samples<2:
        raise Exception(f"condpair has not enough samples: c1:{len(samples_c1)} c2: {len(samples_c2)}, skipping")

    minrep_c1 = get_minrep_for_cond(samples_c1, minrep)
    minrep_c2 = get_minrep_for_cond(samples_c2, minrep)
    df_c1 = unnormed_df.loc[:, samples_c1].dropna(thresh=minrep_c1, axis=0)
    df_c2 = unnormed_df.loc[:, samples_c2].dropna(thresh=minrep_c2, axis=0)
    if (len(df_c1.index)<5) | (len(df_c2.index)<5):
        raise Exception(f"condpair has not enough data for processing c1: {len(df_c1.index)} c2: {len(df_c2.index)}, skipping")

    return df_c1, df_c2

def get_minrep_for_cond(c_samples, minrep):
    if minrep is None: #in the case of None, no nans will be allowed
        return None
    num_samples = len(c_samples)
    if num_samples<minrep:
        return num_samples
    else:
        return minrep