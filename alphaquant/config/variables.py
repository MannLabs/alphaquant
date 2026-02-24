import pandas as pd

QUANT_ID = "quant_id"
MIN_PVAL = 1e-16
PROGRESS_FOLDER = "progress"
PREFER_PRECURSORS_FOR_CLUSTERING = True
PEPTIDE_OUTLIER_FILTERING = True
PTM_FRAGMENT_SELECTION = False
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

def set_quant_id(quant_id):
    global QUANT_ID
    QUANT_ID = quant_id

def set_peptide_outlier_filtering(peptide_outlier_filtering):
    global PEPTIDE_OUTLIER_FILTERING
    PEPTIDE_OUTLIER_FILTERING = peptide_outlier_filtering

def set_ptm_fragment_selection(is_ptm: bool):
    global PTM_FRAGMENT_SELECTION
    PTM_FRAGMENT_SELECTION = bool(is_ptm)

# Backwards-compat alias
def set_phospho_fragment_selection(is_phospho: bool):
    set_ptm_fragment_selection(is_phospho)


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
