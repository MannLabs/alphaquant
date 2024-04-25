import pandas as pd

QUANT_ID = "quant_id"
MIN_PVAL = 1e-16
PROGRESS_FOLDER = "progress"

def determine_variables(input_file):
    _determine_quant_id(input_file)


def _determine_quant_id(input_file):
    global QUANT_ID
    if "aq_reformat.tsv" in input_file:
        input_df = pd.read_csv(input_file, sep="\t", nrows=3)        
        if "quant_id" in input_df.columns:
            QUANT_ID = "quant_id"
        elif "ion" in input_df.columns:
            QUANT_ID = "ion"

def set_quant_id(quant_id):
    global QUANT_ID
    QUANT_ID = quant_id