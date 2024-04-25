import logging
import alphabase.quantification.quant_reader.config_dict_loader as ab_config_loader
import os
import re

############################################
# Set up logging

def setup_logging():
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
setup_logging()

############################################
# Adapt some functions in the alphabase quant reader package

def _get_seperator(input_file):
    filename = str(input_file)
    sep = None
    if '.csv' in filename:
        sep = ','
    elif '.tsv' in filename or '.txt' in filename:
        sep = '\t'

    if sep is None:
        sep = "\t"
    return sep

def _get_original_file_from_aq_reformat(input_file):
    input_file = os.path.abspath(input_file)
    path_elements = input_file.split(os.sep)
    filename = path_elements[-1]
    filename_new = re.match("(.*)(\..*\.)(aq_reformat\.tsv)",filename).group(1)
    path_elements[-1] = filename_new
    path_elements.remove("progress") #remove the progress folder
    return os.sep.join(path_elements)

# Override the function in the package
ab_config_loader._get_seperator = _get_seperator
ab_config_loader._get_original_file_from_aq_reformat = _get_original_file_from_aq_reformat