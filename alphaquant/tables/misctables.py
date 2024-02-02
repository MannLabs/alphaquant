import pandas as pd
import os
import alphabase.quantification.quant_reader.config_dict_loader as abconfigloader
import alphaquant.utils.utils as aqutils

class RunConfigTableCreator():
    def __init__(self, runconfig):
        self._runconfig = runconfig

        self.runconfig_df = None

        self._define_results_df()

    def _define_results_df(self):
        method_params = self._get_methods_dict_from_runconfig()
        self.runconfig_df = pd.Series(method_params)

    def _get_methods_dict_from_runconfig(self):
        method_params = {}
        local_vars = self._runconfig.__dict__
        for x in local_vars.keys():
            if local_vars[x] is None:
                continue
            if isinstance(local_vars[x], pd.DataFrame):
                continue

            if (("_df" not in x) and ('condpair' not in x) and ('sys'!=x) and ('runconfig' != x)):
                if ("input_file" in x) or ("results_dir" in x):
                    method_params[x] = os.path.abspath(local_vars[x])
                else:
                    method_params[x] = local_vars[x]
        return method_params


class AnnotationFileCreator():
    """Creates a file that can be used for annotation of the results table."""
    def __init__(self, input_file, input_type_to_use, annotation_columns_user_defined = None):
        self.annotation_filename = None

        self._input_file = input_file
        self._input_type_to_use = input_type_to_use
        self._annotation_columns_user_defined = annotation_columns_user_defined
        self._id_columns = None
        self._annotation_columns = None
        self._annotation_df = None
        
        self._define_id_and_annotation_columns()
        self._define_and_write_annotation_file_if_applicable()
    
    def _define_id_and_annotation_columns(self):
        config_dict = abconfigloader.get_input_type_and_config_dict(input_file = self._input_file, input_type_to_use = self._input_type_to_use)[1] #load the config dict which contains the relevant columns
        self._id_columns = config_dict.get("protein_cols")
        if self._annotation_columns_user_defined is not None:
            self._annotation_columns = self._annotation_columns_user_defined
        else:
            self._annotation_columns = config_dict.get("annotation_columns", None)

    def _define_and_write_annotation_file_if_applicable(self):
        if self._annotation_columns is not None:
            self._define_annotation_df()
            self._define_annotation_filename()
            self._write_annotation_file()
    
    def _define_annotation_df(self):
        self._annotation_df = pd.read_csv(self._input_file, usecols = self._id_columns + self._annotation_columns, sep = "\t",encoding ='latin1').drop_duplicates()
        self._annotation_df["protein"] = self._annotation_df[self._id_columns].apply(lambda x : "_".join(x), axis = 1)
        self._annotation_df = self._annotation_df.drop(columns = self._id_columns)
    
    def _define_annotation_filename(self):
        filename_shortened = aqutils.remove_file_extension(self._input_file)
        self.annotation_filename = f"{filename_shortened}.annotation.tsv"

    def _write_annotation_file(self):
        self._annotation_df.to_csv(self.annotation_filename, sep = "\t", index = None)
    
        