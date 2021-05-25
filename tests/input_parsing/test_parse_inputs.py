import yaml
import pandas as pd
from matplotlib_venn import venn2

def compare_generic_table_with_original(preprocessed_input_df, original_input_file, config_yaml,input_typename_config, sep = "\t"):
    id2quant_orig, id2quant_preproc = get_processed_original_id2quant_maps(preprocessed_input_df, original_input_file, config_yaml,input_typename_config)
    venn2([set(id2quant_orig.keys()), set(id2quant_preproc.keys())])


def get_processed_original_id2quant_maps(preprocessed_input_df, original_input_file, config_yaml,input_typename_config, sep = "\t"):
    config_all = yaml.safe_load(open(config_yaml, 'r'))
    config_dict = config_all.get(input_typename_config)
    id_cols = config_dict.get("ion_cols") + [config_dict.get("sample_ID")]
    quant_col = config_dict.get("quant_ID")
    id2quant_orig = get_id2quant_original(original_input_file, id_cols, quant_col, sep)
    id2quant_preproc = get_id2quant_processed(preprocessed_input_df, id_cols, quant_col)
    print(id2quant_preproc)
    return id2quant_orig, id2quant_preproc


def get_id2quant_original(original_input_file, id_cols, quant_col, sep):
    print(id_cols)
    print(id_cols+[quant_col])
    orig_df = pd.read_csv(original_input_file, sep=sep, usecols= id_cols+[quant_col])
    orig_df["compareID"] = orig_df[id_cols].sum(axis = 1)
    display(orig_df)
    id2quant = dict(zip(orig_df["compareID"], orig_df[quant_col]))
    return id2quant


def get_id2quant_processed(preprocessed_input_df, id_cols, quant_col):
    compare_IDs = []
    quantvals = []
    for column in preprocessed_input_df.columns:
        if(column == "protein"):
            continue
        ids = pd.Series([column for x in range(len(preprocessed_input_df.index))])
        compare_IDs.extend(list(preprocessed_input_df.index + ids))
        quantvals.extend(preprocessed_input_df[column])
    id2quant = dict(zip(ids, quantvals))
    return id2quant


import alphaquant.diffquant_utils as aqutils
import os
tabledir = os.path.join("..", "..", "test_data", "input_table_formats")
results_folder = os.path.join(".", "results")
print(os.path.abspath(results_folder))

input_file = os.path.join(tabledir,"diann.tsv" )
samplemap_file = os.path.join(tabledir, "samplemap.diann.tsv")

input_data = aqutils.import_data(input_file, results_folder=results_folder)
samplemap_df = aqutils.load_samplemap(samplemap_file)
input_processed, samplemap_df_processed = aqutils.prepare_loaded_tables(input_data, samplemap_df)

compare_generic_table_with_original(input_processed, input_file, os.path.join("..", "..", "intable_config.yaml"), "diann_precursor")