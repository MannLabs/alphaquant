import alphaquant.diffquant.diffutils as aqutils
import alphaquant.config.config as aqconfig
import alphaquant.plotting.classify as aq_plot_classify
import alphaquant.classify.ml_info_table as aq_ml_info_table
import alphaquant.config.variables as aq_conf_vars
import alphaquant.classify.classification_utils as aq_class_utils
import alphaquant.classify.training_functions as aq_class_train

import sklearn
import numpy as np
import anytree
import pandas as pd

import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)


def assign_predictability_scores_stacked(protein_nodes, results_dir, ml_info_file ,name,samples_used, min_num_precursors=3, prot_fc_cutoff =  0, replace_nans = False,
                                         plot_predictor_performance = False, performance_metrics = {}, shorten_features_for_speed = True):
    #protnorm peptides should always be true, except when the dataset run tests different injection amounts

    #add predictability scores to each precursor
    #prepare the input table with all the relevant features for machine learning

    protein_nodes = list(sorted(protein_nodes, key  = lambda x : x.name))
    
    precursor_selector = PrecursorForTrainingSelector(protein_nodes, min_num_precursors = min_num_precursors, prot_fc_cutoff = prot_fc_cutoff)
    
    LOGGER.info(f"{precursor_selector.num_precursors_suitable_for_training} of {precursor_selector.num_precursors_total} selected for training")

    if len(precursor_selector.precursors_suitable_for_training)<100:
        LOGGER.info(f"too few precursors suitable for training, skipping ml")
        return False


    
    acquisition_info_df = aq_ml_info_table.MLInfoTableLoader(ml_info_file, samples_used).ml_info_df

    ml_input_for_training = MLInputTableCreator(precursor_selector.precursors_suitable_for_training, acquisition_info_df, define_y = True, replace_nans = replace_nans)
    ml_input_remaining = MLInputTableCreator(precursor_selector.precursors_not_suitable_for_training, acquisition_info_df, define_y = False, replace_nans = replace_nans)
    align_ml_input_tables_if_necessary(ml_input_for_training, ml_input_remaining)


    featurenames_str = ', '.join(ml_input_for_training.featurenames)
    LOGGER.info(f"starting ML prediction using features {featurenames_str}")
    models, test_set_predictions, y_pred_cv = aq_class_train.train_fast_gradient_boosting(ml_input_for_training.X, ml_input_for_training.y,
                                                                           num_splits=5, shorten_features_for_speed=False)

    # Use out-of-fold predictions for ml_input_for_training.X
    y_pred = y_pred_cv

    y_pred_remaining = predict_on_models(models, ml_input_remaining.X)
    y_pred_total = np.concatenate([y_pred, y_pred_remaining])
    ml_scores = convert_y_pred_to_ml_score(y_pred_total) #convert to 0-1 scale where higher is better

    performance_metrics["r2_score"] = sklearn.metrics.r2_score(ml_input_for_training.y, y_pred)

    LOGGER.info("performed ML prediction")

    #define plot outdir
    results_dir_plots =f"{results_dir}/{name}"
    aqutils.make_dir_w_existcheck(results_dir_plots)
    if plot_predictor_performance:

        aq_plot_classify.scatter_ml_regression_testsets(test_set_predictions, results_dir_plots)
        aq_plot_classify.plot_feature_importance_per_model(models, ml_input_for_training.featurenames, 10, results_dir_plots)
        aq_plot_classify.plot_value_histogram(ml_scores, results_dir_plots)



    
    ionnames_total = ml_input_for_training.ionnames + ml_input_remaining.ionnames
    all_precursors = precursor_selector.precursors_suitable_for_training + precursor_selector.precursors_not_suitable_for_training
    
    #annotate the precursor nodes
    aq_class_utils.annotate_precursor_nodes( ml_scores, ionnames_total, all_precursors) #two new variables added to each node:
    return True



class PrecursorForTrainingSelector:
    def __init__(self, protein_nodes,  min_num_precursors=3, prot_fc_cutoff =  0.75):

        self._protein_nodes = protein_nodes
        self._precursor_cutoff = min_num_precursors
        self._prot_fc_cutoff = prot_fc_cutoff

        self.precursors_suitable_for_training = [] # the precursors that are used for training the ML model
        self.precursors_not_suitable_for_training = [] # the precursors that are not used for training the ML model

        self._select_precursors_for_training()

        self.num_precursors_suitable_for_training = len(self.precursors_suitable_for_training)
        self.num_precursors_not_suitable_for_training = len(self.precursors_not_suitable_for_training)
        self.num_precursors_total = self.num_precursors_suitable_for_training + self.num_precursors_not_suitable_for_training

    
    def _select_precursors_for_training(self):
        for protein_node in self._protein_nodes:  
            precursors = self._get_precursors(protein_node)              
            if (abs(protein_node.fc) < self._prot_fc_cutoff) or (len(precursors)<self._precursor_cutoff):
                self.precursors_not_suitable_for_training.extend(precursors)
            else:
                self.precursors_suitable_for_training.extend(precursors)
    
    @staticmethod
    def _get_precursors(protein_node):
        if protein_node.leaves[0].parent.level == "mod_seq":
            return protein_node.leaves
        else:
            return anytree.findall(protein_node, filter_= lambda x : (x.level == "mod_seq_charge"))


class MLInputTableCreator:
    def __init__(self, precursors, acquisition_info_df, define_y,replace_nans = False, numeric_threshold = 0.99):
        self._precursors = precursors
        self._acquisition_info_df = acquisition_info_df
        self._replace_nans = replace_nans
        self._numeric_threshold = numeric_threshold

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
        if define_y: # y should only be defined if the protein has enough precursors in order to get a meaningful y
            self._define_y()

    def _define_merged_df(self):
        node_features_df = aq_class_utils.collect_node_parameters(self._precursors)
        self._merged_df = aqutils.merge_acquisition_df_parameter_df(self._acquisition_info_df, node_features_df)

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
        ion2fc = {x.name: self._get_protnormed_fc(x) for x in self._precursors}
        self.y = np.array([ion2fc.get(ion) for ion in self.ionnames])
    
    @staticmethod
    def _get_protnormed_fc(precursor_node):
        fc_precursor = precursor_node.fc
        protein_node = aq_class_utils.find_node_parent_at_level(precursor_node, "gene")
        protein_fc = protein_node.fc
        return fc_precursor - protein_fc


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



def predict_on_models(models, X):
    y_preds = [model.predict(X) for model in models]
    y_pred = np.mean(y_preds, axis=0)
    return y_pred

def convert_y_pred_to_ml_score(y_pred):
    """
    Convert model predictions to normalized ML scores.

    This function takes the absolute values of the predictions, normalizes them
    to a 0-1 scale, and then inverts the scale so that higher scores represent
    better predictions.

    Args:
        y_pred (numpy.ndarray): An array of model predictions. Can contain
            positive or negative float values.

    Returns:
        numpy.ndarray: An array of ML scores normalized between 0 and 1, where
            higher scores indicate better predictions. If all input predictions
            are 0, the function returns an array of 1s.

    Note:
        - If the input array contains all zeros, the output will be an array of ones.
        - The function assumes that lower original predictions are better, and
          the output inverts this so that higher ML scores are better.
    """

    # Take the absolute values to ensure all predictions are non-negative
    abs_values = np.abs(y_pred)

    # Normalize these values to a 0-1 scale
    max_value = np.max(abs_values)
    if max_value > 0:
        normalized_scores = abs_values / max_value
    else:
        normalized_scores = abs_values  # If max_value is 0, all values are zero and hence already normalized

    # Reverse the order so that higher scores are better
    ml_score = 1 - normalized_scores

    return ml_score




