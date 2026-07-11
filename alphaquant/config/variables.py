import pandas as pd

QUANT_ID = "quant_id"
MIN_PVAL = 1e-16
PROGRESS_FOLDER = "progress"
PREFER_PRECURSORS_FOR_CLUSTERING = True
PEPTIDE_OUTLIER_FILTERING = True
PEPTIDE_OUTLIER_REGULATION_NORMALIZATION_FACTOR = 15.0
OUTLIER_CORRECTION_FACTOR = 1.0
PTM_FRAGMENT_SELECTION = False
MAX_N_FRAGMENTS = None
ION_OUTLIER_MAD_THRESHOLD = None
CLASSIC_FRAGMENT_OUTLIER_FILTERING = False
ICC_NULL_PVAL_THRESHOLD = 0.1
NUM_BG_CONTEXTS = 10
# When residual-decorrelation pruning collapses a parent to a single surviving
# child, aggregate that parent via the MEDIAN of ALL its children's z-values
# (shared-signal estimate for near-duplicate siblings) instead of using the lone
# survivor. Applied at every tree level. Default off. See set_median_on_collapse.
MEDIAN_ON_COLLAPSE = False
# When True, residual-decorrelation records the mean pairwise correlation among the
# SURVIVING children of each parent (the residual, post-pruning correlation) into
# node.icc_correction, so the Stouffer aggregation applies a design-effect
# deff=1+(n-1)*rho. This corrects the between-peptide correlation that pruning cannot
# remove (homogeneous, no droppable subset). Where pruning already decorrelated the
# children the residual rho ~ 0, so the correction is a no-op. Default off.
RESIDUAL_DEFF_CORRECTION = False
# Maximum number of peptides combined per protein (closest to median z are kept).
# Bounds significance for very peptide-rich proteins. Set to None to disable the cap.
MAX_PEPTIDES_PER_PROTEIN = 31
CONDITION_PAIR_SEPARATOR = "_VS_"

#prefixes for the different ion types
SEQ = "SEQ"
MOD = "MOD"
CHARGE = "CHARGE"
FRG = "FRG"
ION = "ION"

INPUT_TYPE = None   # e.g. "diann_precursor_fragion", set via set_input_config()
CONFIG_DICT = None  # the full config dict for the detected input type


def determine_variables(input_file, input_type):
    _determine_quant_id(input_file)
    _determine_prefer_precursors_for_clustering(input_type)


def _determine_quant_id(input_file):
    global QUANT_ID
    if "aq_reformat.tsv" in input_file:
        input_df = pd.read_csv(input_file, sep="\t", nrows=3)
        if "quant_id" in input_df.columns:
            QUANT_ID = "quant_id"
        elif "ion" in input_df.columns:
            QUANT_ID = "ion"

def _determine_prefer_precursors_for_clustering(input_type):
    global PREFER_PRECURSORS_FOR_CLUSTERING
    if "precursor_fragion" in input_type:
        PREFER_PRECURSORS_FOR_CLUSTERING = True
    else:
        PREFER_PRECURSORS_FOR_CLUSTERING = False

def set_peptide_outlier_filtering(peptide_outlier_filtering):
    global PEPTIDE_OUTLIER_FILTERING
    PEPTIDE_OUTLIER_FILTERING = peptide_outlier_filtering

def set_median_on_collapse(median_on_collapse):
    global MEDIAN_ON_COLLAPSE
    MEDIAN_ON_COLLAPSE = bool(median_on_collapse)

def set_residual_deff_correction(residual_deff_correction):
    global RESIDUAL_DEFF_CORRECTION
    RESIDUAL_DEFF_CORRECTION = bool(residual_deff_correction)

def set_max_peptides_per_protein(max_peptides_per_protein):
    global MAX_PEPTIDES_PER_PROTEIN
    MAX_PEPTIDES_PER_PROTEIN = (int(max_peptides_per_protein)
                                if max_peptides_per_protein is not None else None)

def set_outlier_correction_factor(outlier_correction_factor):
    global OUTLIER_CORRECTION_FACTOR
    OUTLIER_CORRECTION_FACTOR = float(outlier_correction_factor)

def set_max_n_fragments(max_n_fragments):
    global MAX_N_FRAGMENTS
    MAX_N_FRAGMENTS = int(max_n_fragments) if max_n_fragments is not None else None

def set_ion_outlier_mad_threshold(threshold):
    global ION_OUTLIER_MAD_THRESHOLD
    ION_OUTLIER_MAD_THRESHOLD = float(threshold) if threshold is not None else None

def set_classic_fragment_outlier_filtering(enabled):
    global CLASSIC_FRAGMENT_OUTLIER_FILTERING
    CLASSIC_FRAGMENT_OUTLIER_FILTERING = bool(enabled)

def set_icc_null_pval_threshold(threshold):
    """Set the ICC null p-value threshold.

    Args:
        threshold: Float applied to all ICC-estimation levels.
    """
    global ICC_NULL_PVAL_THRESHOLD
    ICC_NULL_PVAL_THRESHOLD = float(threshold)


def get_icc_null_pval_threshold(node_type=None):
    """Return the ICC null p-value threshold for *node_type*.

    ``node_type`` is accepted for compatibility with callers that ask for the
    level-specific value, but one global threshold is now used for all levels.
    """
    return ICC_NULL_PVAL_THRESHOLD

def set_ptm_fragment_selection(is_ptm: bool):
    global PTM_FRAGMENT_SELECTION
    PTM_FRAGMENT_SELECTION = bool(is_ptm)

def set_input_config(input_type, config_dict):
    """Store the detected input type and its full config dict as module globals.

    Called once during pipeline setup so that other modules (e.g.
    ``background_distributions``) can inspect the active configuration
    without passing it through every function signature.

    Args:
        input_type (str): Identifier of the detected input format
            (e.g. ``"diann_precursor_fragion"``).
        config_dict (dict): The complete YAML config dict for *input_type*.

    Side effects:
        Sets ``INPUT_TYPE`` and ``CONFIG_DICT`` at module level.
    """
    global INPUT_TYPE, CONFIG_DICT
    INPUT_TYPE = input_type
    CONFIG_DICT = config_dict
