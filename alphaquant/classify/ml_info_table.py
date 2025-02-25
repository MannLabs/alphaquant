import alphabase.quantification.quant_reader.config_dict_loader as ab_config_loader
import alphabase.quantification.quant_reader.table_reformatter as ab_table_reformatter
import alphaquant.utils.reader_utils as aq_reader_utils

import pandas as pd
import alphaquant.utils.utils as aq_utils
import numpy as np
import itertools

class MLInfoTableCreator():
    def __init__(self, input_file, input_type_to_use, modfication_type):
        self._input_file = input_file
        self._input_type_to_use = input_type_to_use
        self._modification_type = modfication_type

        self._ml_info_df = None

        self.ml_info_filename = None

        self._define_ml_info_df()
        self._define_ml_info_filename()
        self._write_ml_info_table()

    def _define_ml_info_df(self):

        """We need to stitch together the relevant columns to an ionname. For this we use the 'config dict' from the quant table reader. An example config_dict is given below

            ion_hierarchy:
                fragion:
                order: [SEQ, MOD, CHARGE, FRGION]
                mapping:
                    SEQ:
                    - Stripped.Sequence
                    MOD:
                    - Modified.Sequence
                    CHARGE:
                    - Precursor.Charge
                    FRGION:
                    - Fragment.Quant.Corrected
        """
        input_df = aq_reader_utils.read_file(self._input_file, sep="\t")
        _, config_dict, _ = ab_config_loader.get_input_type_and_config_dict(self._input_file, self._input_type_to_use)
        ion_hierarchy = config_dict.get("ion_hierarchy")
        #get first dict entry
        hierarchy_dict = next(iter(ion_hierarchy.values())) #just get the first entry, as we only want SEQ, MOD, CHARGE, which should be the same in all cases
        ion_order = ["SEQ", "MOD", "CHARGE"]
        ion_mapping = hierarchy_dict.get("mapping")
        ion_headers_merged, ion_headers_grouped = ab_table_reformatter.get_ionname_columns(ion_mapping, ion_order)
        ion_headers_grouped = [ion_mapping.get("SEQ"), ion_mapping.get("MOD"), ion_mapping.get("CHARGE")]
        ion_headers_merged = list(itertools.chain(*ion_headers_grouped))

        input_df = ab_table_reformatter.add_index_and_metadata_columns(df_subset=input_df, ion_hierarchy_local=ion_order, ion_headers_grouped=ion_headers_grouped, quant_id_dict=None, hierarchy_type=None)
        input_df = input_df.reset_index().drop(columns=ion_headers_merged)
        sample_column = config_dict.get("sample_ID")
        input_df = input_df.rename(columns={sample_column: "sample_ID"})
        input_df = input_df.set_index(["quant_id", "sample_ID"])
        self._ml_info_df = input_df.select_dtypes(include=[np.number])
        self._ml_info_df = self._ml_info_df.reset_index()
        if self._modification_type is not None:
            self._adapt_precursor_name_to_modification_type()

    def _adapt_precursor_name_to_modification_type(self):
        self._ml_info_df["quant_id"] = self._ml_info_df["quant_id"].apply(lambda x: x.replace(self._modification_type, ""))


    def _define_ml_info_filename(self):
        self.ml_info_filename = aq_utils.get_progress_folder_filename(self._input_file, ".ml_info_table.tsv")

    def _write_ml_info_table(self):
        self._ml_info_df.to_csv(self.ml_info_filename, sep="\t", index=False)


class MLInfoTableLoader():
    def __init__(self, ml_info_file, samples_used):
        self.ml_info_df = aq_reader_utils.read_file(ml_info_file, sep="\t")

        self._samples_used = samples_used

        self._subset_df_to_relevant_samples()

    def _subset_df_to_relevant_samples(self):
        self.ml_info_df = self.ml_info_df[self.ml_info_df["sample_ID"].isin(self._samples_used)]
        self.ml_info_df = self.ml_info_df.drop(columns=["sample_ID"])

