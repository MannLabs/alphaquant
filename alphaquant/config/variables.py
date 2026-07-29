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

MEDIAN_ON_COLLAPSE = False


RESIDUAL_DEFF_CORRECTION = True

RESIDUAL_DEFF_SMALLN_TOTAL = 7

RESIDUAL_DECORR_CORR_MODE = "cap"
RESIDUAL_DECORR_CORR_CAP = 10

MAX_PEPTIDES_PER_PROTEIN = None
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


def set_residual_deff_smalln_total(residual_deff_smalln_total):
    global RESIDUAL_DEFF_SMALLN_TOTAL
    RESIDUAL_DEFF_SMALLN_TOTAL = int(residual_deff_smalln_total) if residual_deff_smalln_total else 0


def set_residual_decorr_corr_mode(residual_decorr_corr_mode):
    global RESIDUAL_DECORR_CORR_MODE
    RESIDUAL_DECORR_CORR_MODE = str(residual_decorr_corr_mode) if residual_decorr_corr_mode else "cap"


def set_residual_decorr_corr_cap(residual_decorr_corr_cap):
    global RESIDUAL_DECORR_CORR_CAP
    RESIDUAL_DECORR_CORR_CAP = int(residual_decorr_corr_cap) if residual_decorr_corr_cap else 10

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
