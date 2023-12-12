import pandas as pd


class MedianConditionManager():
    def __init__(self, input_file, samplemap_file):

        self._input_file = input_file
        self._samplemap_file = samplemap_file
        self._samplemap_filename_adapted = None

        self.input_filename_adapted = None
        self._median_condition_creator = None
        self.samplemap_df_adapted = None

        self._define_median_condition_creator()
        self._define_adapted_filenames()
        self._save_adapted_files()
        self._define_adapted_samplemap_df()
    
    def _define_median_condition_creator(self):

        input_df = pd.read_csv(self._input_file, sep = "\t")
        samplemap_df = pd.read_csv(self._samplemap_file, sep = "\t")
        self._median_condition_creator = MedianConditionCreator(samplemap_df, input_df)
    
    def _define_adapted_filenames(self):
        self._samplemap_filename_adapted = self._adapt_filename(self._samplemap_file)
        self.input_filename_adapted = self._adapt_filename(self._input_file)

    @staticmethod
    def _adapt_filename(filename):
        return filename.replace(".tsv", "_w_median.tsv")
    
    def _save_adapted_files(self):
        self._median_condition_creator.extended_input_df.to_csv(self.input_filename_adapted, sep = "\t")
        self._median_condition_creator.extended_samplemap_df.to_csv(self._samplemap_filename_adapted, sep = "\t", index = None)

    def _define_adapted_samplemap_df(self):
        self.samplemap_df_adapted = self._median_condition_creator.extended_samplemap_df




class MedianConditionCreator():
    def __init__(self, samplemap_df, input_df_aqformat):
        self._samplemap_df = samplemap_df
        self._input_df_aqformat = input_df_aqformat.set_index(["protein", "quant_id"])

        self._number_replicates = self._determine_number_replicates()

        self.extended_input_df = self._define_extended_input_df()
        self.extended_samplemap_df = self._define_extended_samplemap_df()

    def _define_extended_input_df(self):
        input_df_subset = self._subset_input_df_to_samplemap()
        median_sample_df = self._define_median_sample_df()
        return pd.concat([input_df_subset, median_sample_df], axis="columns")

    def _subset_input_df_to_samplemap(self):
        return self._input_df_aqformat[self._samplemap_df["sample"].to_list()]
    
    def _define_median_sample_df(self):
        replicate_intensities = []
        for idx in range(self._number_replicates):
            replicate_intensities.append(self._get_median_vals_for_replicate_idx(idx))
        median_sample_df = pd.concat(replicate_intensities, axis="columns")
        median_sample_df.columns = [f"median_rep{idx}" for idx in range(self._number_replicates)]
        return median_sample_df
    
    def _determine_number_replicates(self):
        replicate_numbers = self._samplemap_df.groupby("condition").size()
        return replicate_numbers.min()
    
    def _get_median_vals_for_replicate_idx(self, replicate_idx):
        list_of_sample_intensities = []
        for (groupname, cond_df) in self._samplemap_df.groupby("condition"):
            expname = cond_df["sample"].to_list()[replicate_idx]
            sample_intensities = self._input_df_aqformat[expname]
            list_of_sample_intensities.append(sample_intensities)
        selected_intensities_df =  pd.concat(list_of_sample_intensities, axis="columns")
        median_intensities = selected_intensities_df.median(axis="columns")

        return median_intensities
    
    def _define_extended_samplemap_df(self):
        df_to_extend = pd.DataFrame({'sample' : [f"median_rep{idx}" for idx in range(self._number_replicates)] , 'condition' : ['median_reference' for idx in range(self._number_replicates)]})
        return pd.concat([self._samplemap_df, df_to_extend], axis="rows")
    

def get_all_conds_relative_to_median(samplemap_df):
    conds = samplemap_df["condition"].unique()
    condpair_combinations = [(x, "median_reference") for x in conds]
    return condpair_combinations
