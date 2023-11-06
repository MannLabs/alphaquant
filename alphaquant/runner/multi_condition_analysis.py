import pandas as pd



def add_and_save_median_condition(input_file, samplemap_file):
    input_df = pd.read_csv(input_file, sep = "\t")
    samplemap_df = pd.read_csv(samplemap_file, sep = "\t")

    median_condition_creator = MedianConditionCreator(samplemap_df, input_df)

    median_condition_creator.extended_input_df.to_csv(input_file, sep = "\t")
    median_condition_creator.extended_samplemap_df.to_csv(samplemap_file, sep = "\t", index = None)




class MedianConditionCreator():
    def __init__(self, samplemap_df, input_df_aqformat):
        self._samplemap_df = samplemap_df
        self._input_df_aqformat = input_df_aqformat.set_index(["protein", "quant_id"])

        self._number_replicates = self._determine_number_replicates()

        self.extended_input_df = self._define_extended_input_df()
        self.extended_samplemap_df = self._define_extended_samplemap_df()

    def _define_extended_input_df(self):
        median_sample_df = self._define_median_sample_df()
        return pd.concat([self._input_df_aqformat, median_sample_df], axis="columns")

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