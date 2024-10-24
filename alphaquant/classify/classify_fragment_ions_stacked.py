import alphaquant.diffquant.diffutils as aqutils
import alphaquant.classify.classify_ions as aq_class_ions
import alphaquant.config.config as aqconfig
import alphaquant.plotting.classify as aq_plot_classify

import numpy as np
import anytree
import alphaquant.config.variables as aq_conf_vars
import pandas as pd
import alphaquant.classify.classification_utils as aq_class_utils
import alphaquant.cluster.cluster_utils as aq_cluster_utils

import sklearn.ensemble
import sklearn.linear_model
import sklearn.impute
import sklearn.metrics
import sklearn.model_selection




import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)


def assign_predictability_scores_stacked(protein_nodes, results_dir, name, acquisition_info_df, min_num_fragions=3,  replace_nans = False,
                                         plot_predictor_performance = False, performance_metrics = {}, shorten_features_for_speed = False):
    #protnorm peptides should always be true, except when the dataset run tests different injection amounts

    #add predictability scores to each fragion
    #prepare the input table with all the relevant features for machine learning

    protein_nodes = list(sorted(protein_nodes, key  = lambda x : x.name))
    
    fragion_selector = FragionForTrainingSelector(protein_nodes, min_num_fragions = min_num_fragions)
    
    LOGGER.info(f"{fragion_selector.num_fragions_suitable_for_training} of {fragion_selector.num_fragions_total} selected for training")

    if fragion_selector.num_fragions_suitable_for_training<100:
        LOGGER.info(f"too few fragions suitable for training, skipping ml")
        return False




    ml_input_for_training = MLInputTableCreatorFragions(fragion_selector.fragions_suitable_for_training, acquisition_info_df, define_y = True, replace_nans = replace_nans)
    ml_input_remaining = MLInputTableCreatorFragions(fragion_selector.fragions_not_suitable_for_training, acquisition_info_df, define_y = False, replace_nans = replace_nans)
    align_ml_input_tables_if_necessary(ml_input_for_training, ml_input_remaining)


    featurenames_str = ', '.join(ml_input_for_training.featurenames)
    LOGGER.info(f"starting RF prediction using features {featurenames_str}")
    
    models = train_random_forest_ensemble(ml_input_for_training.X, ml_input_for_training.y, num_splits = 5, shorten_features_for_speed=shorten_features_for_speed)

    y_pred = predict_on_models(models,ml_input_for_training.X)
    y_pred_remaining = predict_on_models(models, ml_input_remaining.X)
    y_pred_total = np.concatenate([y_pred, y_pred_remaining])

    performance_metrics["r2_score"] = sklearn.metrics.r2_score(ml_input_for_training.y, y_pred)

    LOGGER.info("performed RF prediction")

    #define plot outdir
    results_dir_plots =f"{results_dir}/{name}"
    aqutils.make_dir_w_existcheck(results_dir_plots)
    if plot_predictor_performance:

        aq_plot_classify.scatter_ml_regression(ml_input_for_training.y, y_pred, results_dir_plots)
        #aq_plot_classify.compute_and_plot_feature_importances_stacked_rf(model=stacked_regressor, X_val=ml_input_for_training.X, y_val=ml_input_for_training.y, feature_names=ml_input_for_training.featurenames, top_n=10, results_dir=results_dir_plots)
        aq_plot_classify.plot_feature_importance_per_model(models, ml_input_for_training.featurenames, 10, results_dir_plots)
        aq_plot_classify.plot_value_histogram(y_pred_total, results_dir_plots)


    
    ionnames_total = ml_input_for_training.ionnames + ml_input_remaining.ionnames
    all_fragion_basenodes = fragion_selector.fragions_suitable_for_training + fragion_selector.fragions_not_suitable_for_training
    
    #annotate the fragion nodes
    annotate_fragion_basenodes(all_fragion_basenodes, ionnames_total, y_pred_total) #two new variables added to each node:
    #update_fold_change_of_the_fragion_iontype_node(fragion_selector.fragment_iontype_nodes)
    propagate_new_fcs_along_the_tree(protein_nodes)

    # re_order_fragion_iontype_nodes_by_score(fragion_selector.fragment_iontype_nodes)
    # propagate_new_clusters_along_the_tree(protein_nodes)

    return True



class FragionForTrainingSelector:
    def __init__(self, protein_nodes,  min_num_fragions = 3 ):

        self._protein_nodes = protein_nodes
        self._min_num_fragions = min_num_fragions
        self.fragment_iontype_nodes = []
        self.non_fragment_iontype_nodes = []

        self.fragions_suitable_for_training = [] # the fragions that are used for training the ML model
        self.fragions_not_suitable_for_training = [] # the fragions that are not used for training the ML model

        self._define_iontype_nodes()
        self._select_fragions_for_training()

        self.num_fragions_suitable_for_training = len(self.fragions_suitable_for_training)
        self.num_fragions_not_suitable_for_training = len(self.fragions_not_suitable_for_training)
        self.num_fragions_total = self.num_fragions_suitable_for_training + self.num_fragions_not_suitable_for_training

    
    def _define_iontype_nodes(self):
        for protein_node in self._protein_nodes:
            iontype_nodes = anytree.search.findall(protein_node, filter_=lambda node: node.level == "ion_type")
            for iontype_node in iontype_nodes:
                if iontype_node.type == "frgion":
                    self.fragment_iontype_nodes.append(iontype_node)
                else:
                    self.non_fragment_iontype_nodes.append(iontype_node)


    def _select_fragions_for_training(self):
        for fragment_node in self.fragment_iontype_nodes:
            if  (len(fragment_node.children)<self._min_num_fragions):
                self.fragions_not_suitable_for_training.extend(fragment_node.children)
            else:
                self.fragions_suitable_for_training.extend(fragment_node.children)


class MLInputTableCreatorFragions:
    def __init__(self, fragions, acquisition_info_df, define_y,replace_nans = False, numeric_threshold = 0.3):
        self._fragions = fragions
        self._acquisition_info_df = acquisition_info_df
        self._replace_nans = replace_nans
        self._numeric_threshold = numeric_threshold #fraction of non-nan values in a column, if less, the column is removed

        self._merged_df = None
        
        self.X = None # the input for the ML model which has corresponding y values, so it is possible to train with this table
        self.y = None
        self.featurenames = None
        self.ionnames = None

        self._define_merged_df()
        self._define_ionnames()
        self._remove_non_numeric_columns_from_merged_df()
        self._define_featurenames()
        self._define_X()
        if define_y: # y should only be defined if the protein has enough fragions in order to get a meaningful y
            self._define_y()

    def _define_merged_df(self):
        node_features_df = self._collect_node_parameters()
        if self._acquisition_info_df is None:
            self._merged_df = node_features_df
            return
        self._merged_df = aqutils.merge_acquisition_df_parameter_df(self._acquisition_info_df, node_features_df)

    def _collect_node_parameters(self):
        node_features = []
        for fragion in self._fragions:
            features_dict = {aq_conf_vars.QUANT_ID : fragion.name,"cv" : fragion.cv, "min_intensity" : fragion.min_intensity, "total_intensity" : fragion.total_intensity, "min_reps" : fragion.min_reps}
            node_features.append(features_dict)
        node_features_df = pd.DataFrame(node_features)
        return node_features_df

    def _define_ionnames(self):
        self.ionnames = list(self._merged_df[aq_conf_vars.QUANT_ID])
        

    def _remove_non_numeric_columns_from_merged_df(self):
        columns_to_drop = []
        self._merged_df = self._merged_df.drop(columns=[aq_conf_vars.QUANT_ID])
        self._merged_df = self._merged_df.apply(lambda col: pd.to_numeric(col, errors='coerce')) #'coerce' will turn non-numeric values into NaN
        
        for column in self._merged_df.columns:
            proportion_non_nans = self._merged_df[column].notna().mean()
            if proportion_non_nans < self._numeric_threshold:
                columns_to_drop.append(column)
        
        self._merged_df = self._merged_df.drop(columns=columns_to_drop)

    def _define_featurenames(self):
        self.featurenames = list(self._merged_df.columns)

    def _define_X(self):
        X_df = self._merged_df.convert_dtypes(convert_integer=True, convert_floating=True)
        if self._replace_nans:
            imputer = sklearn.impute.SimpleImputer(strategy='most_frequent') #stragies are mean, median, most_frequent, constant
            X_imputed = imputer.fit_transform(X_df)
            self.X =  X_imputed
        else:
            self.X = X_df.to_numpy()
    

    def _define_y(self):
        ion2fc = {x.name: self._get_fragnormed_fc(x) for x in self._fragions}
        self.y = np.array([ion2fc.get(ion) for ion in self.ionnames])
    
    @staticmethod
    def _get_fragnormed_fc(base_node):
        base_fc = base_node.fc
        frgion_node = aq_class_utils.find_node_parent_at_type(base_node, "frgion") #the nodes we are classifying are "base" nodes. The "frgion" node has all the base nodes as children
        fragion_fc = frgion_node.fc
        return base_fc - fragion_fc


def align_ml_input_tables_if_necessary(ml_input_1, ml_input_2):
    featurenames_1 = ml_input_1.featurenames
    featurenames_2 = ml_input_2.featurenames

    featurenames_common = list(set(featurenames_1) & set(featurenames_2))
    featurenames_common_ordered = [fn for fn in featurenames_1 if fn in featurenames_common]

    idxs_to_remove_1 = [i for i, featurename in enumerate(featurenames_1) if featurename not in featurenames_common]
    idxs_to_remove_2 = [i for i, featurename in enumerate(featurenames_2) if featurename not in featurenames_common]

    ml_input_1.X = np.delete(ml_input_1.X, idxs_to_remove_1, axis=1)
    ml_input_2.X = np.delete(ml_input_2.X, idxs_to_remove_2, axis=1)

    ml_input_1.featurenames = featurenames_common_ordered
    ml_input_2.featurenames = featurenames_common_ordered


def train_random_forest_ensemble(X, y, shorten_features_for_speed, num_splits=5):
    kf = sklearn.model_selection.KFold(n_splits=num_splits, shuffle=True, random_state=42)
    models = []

    if shorten_features_for_speed:
        max_features = 'sqrt'
    else:
        max_features = 'auto'

    for train_index, _ in kf.split(X):
        X_train, y_train = X[train_index], y[train_index]

        model = sklearn.ensemble.RandomForestRegressor(n_estimators=50,  # Reduced number of trees
                                                       random_state=42,
                                                       n_jobs=-1,  # Utilize all CPU cores
                                                       max_features=max_features)  # Reduce the number of features
        model.fit(X_train, y_train)
        models.append(model)
    
    return models


def predict_on_models(models, X):
    y_preds = [model.predict(X) for model in models]
    y_pred = np.mean(y_preds, axis=0)
    return y_pred



def annotate_fragion_basenodes(all_fragions_basenodes, ionnames_total, y_pred_total):
    ion2pred = dict(zip(ionnames_total, y_pred_total))
    for fragion in all_fragions_basenodes:
        y_pred = ion2pred.get(fragion.name)
        fragion.ml_score_fragion = abs(y_pred)
    


def update_fold_change_of_the_fragion_iontype_node(all_fragions_iontype_nodes): #the iontype nodes are the parents of the basenodes
    for fragion_iontype in all_fragions_iontype_nodes:
        weigths = [2**-abs(fragion.ml_score_fragion) for fragion in fragion_iontype.children]
        fcs = [fragion.fc for fragion in fragion_iontype.children]
        fragion_iontype.fc = np.average(fcs, weights=weigths)
        

def update_fold_change_of_the_mod_seq_ch_node(all_fragions_iontype_nodes):
    for fragion in all_fragions_iontype_nodes:
        fragion.parent.fc = fragion.fc


def propagate_new_fcs_along_the_tree(protein_nodes):
    for protein in protein_nodes:
        for level_nodes in aq_cluster_utils.iterate_through_tree_levels_bottom_to_top(protein):
            nodelevel = level_nodes[0].level
            if nodelevel == "base":
                continue
            if nodelevel == "ion_type":
                continue
            for level_node in level_nodes:
                aq_cluster_utils.aggregate_node_properties(level_node, only_use_mainclust=True, use_fewpeps_per_protein=True)


# def update_nodes_w_ml_score(protnodes):
#     typefilter = globally_initialized_typefilter
#     for prot in protnodes:
#         re_order_depending_on_ml_score(prot, typefilter)


# def re_order_depending_on_ml_score(protnode, typefilter):
#     for idx in range(len(typefilter.type)):
#         type_nodes = anytree.search.findall(protnode, filter_=lambda node: node.type == typefilter.type[idx])
#         if len(type_nodes)==0:
#             continue
#         for type_node in type_nodes: #go through the nodes, re-order the children. Propagate the values from the newly ordered children to the type node
#             child_nodes = type_node.children
#             had_ml_score = hasattr(child_nodes[0], 'ml_score')
#             if had_ml_score:
#                 re_order_clusters_by_ml_score(child_nodes)
#                 aqcluster_utils.aggregate_node_properties(type_node,only_use_mainclust=True, use_fewpeps_per_protein=True)
import copy

def re_order_fragion_iontype_nodes_by_score(fragion_iontype_nodes):
    for fragion_iontype_node in fragion_iontype_nodes:
        re_order_children_of_fragion_iontype_nodes_by_score(fragion_iontype_node.children)


def re_order_children_of_fragion_iontype_nodes_by_score(fragion_basenodes):
    cluster2score = create_cluster2score_dict_fragions(fragion_basenodes)
    clusters = list(cluster2score.keys())
    clusters.sort(key = lambda x : cluster2score.get(x))
    clust2newclust = { clusters[x] :x for x in range(len(clusters))}
    for node in fragion_basenodes:
        cluster_before = copy.copy(node.cluster)
        node.cluster =clust2newclust.get(node.cluster)
        cluster_after = node.cluster
        if cluster_before != cluster_after:
            print(f"cluster changed from {cluster_before} to {cluster_after} for node {node.name}")

def create_cluster2score_dict_fragions(fragion_basenodes): #lower score means better
    cluster2allscores = create_cluster2allscores_dict_fragions(fragion_basenodes)
    cluster2score = {cluster: np.mean(cluster2allscores[cluster]) for cluster in cluster2allscores}
    for cluster, allscores in cluster2allscores.items():
        if len(allscores) < 1:
            cluster2score[cluster] = 10000
    return cluster2score

def create_cluster2allscores_dict_fragions(fragion_basenodes): #lower score means better
    cluster2scores = {}
    for node in fragion_basenodes:
        cluster2scores[node.cluster] = cluster2scores.get(node.cluster, [])
        cluster2scores[node.cluster].append(abs(node.ml_score_fragion))
    return cluster2scores



def propagate_new_clusters_along_the_tree(protein_nodes):
    for protein in protein_nodes:
        for level_nodes in aq_cluster_utils.iterate_through_tree_levels_bottom_to_top(protein):
            nodelevel = level_nodes[0].level
            if nodelevel == "base":
                continue
            for level_node in level_nodes:
                aq_cluster_utils.aggregate_node_properties(level_node, only_use_mainclust=True, use_fewpeps_per_protein=True)


