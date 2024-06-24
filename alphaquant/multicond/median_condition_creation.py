import pandas as pd
import numpy as np


class MedianConditionManager():
    def __init__(self, input_file :str, samplemap_file :str, fraction_missing_values = 0.1):
        """This class reads in a quantitative alphaquant input table and a samplemap file. It calculates new "median intensity" samples,
        adds those to the input table, extends the samplemap accordingly and saves the adapted files under new filenames.
        

        Args:
            input_file (str): alphaquant input file ((protein,quant_id) x sample table in matrix format)
            samplemap_file (str): sample to condition mapping table file
            fraction_missing_values (str): the maximum allowed fraction of missing values allowed when calcuating the median intensity. Defaults to 0.1.
        """

        self._input_df = pd.read_csv(input_file, sep = "\t")
        self._samplemap_df = pd.read_csv(samplemap_file, sep = "\t")
        self._fraction_missing_values = fraction_missing_values
        

        self.input_df_extended = MedianConditionCreator(self._input_df, self._samplemap_df,  self._fraction_missing_values).extended_input_df
        self.samplemap_df_extended = ExtendedSampleMapCreator(self._samplemap_df, self.input_df_extended).extended_samplemap_df
        self.samplemap_filename_adapted = samplemap_file.replace(".tsv", "_w_median.tsv")
        self.input_filename_adapted = input_file.replace(".tsv", "_w_median.tsv")

        self._save_adapted_files()
    
    def _save_adapted_files(self):
        self.samplemap_df_extended.to_csv(self.samplemap_filename_adapted, sep = "\t", index = None)
        self.input_df_extended.to_csv(self.input_filename_adapted, sep = "\t", index = None)



class MedianConditionCreator():
    def __init__(self, input_df_aqformat : pd.DataFrame, samplemap_df : pd.DataFrame, fraction_missing_values : float = 0.1):
        """This class coordinates the calculation of median intensities

        Args:
            input_df_aqformat (pd.DataFrame): (protein,quant_id) x sample dataframe, contaning the intensities
            samplemap_df (pd.DataFrame): sample to condition mapping dataframe
            fraction_missing_values (float, optional): the maximum allowed fraction of missing values allowed when calcuating the median intensity. Defaults to 0.1.
        """
        self._samplemap_df = samplemap_df
        self._input_df_aqformat = np.log2(input_df_aqformat.set_index(["protein", "quant_id"]).replace(0, np.nan))
        self._fraction_missing_values = fraction_missing_values
        self._input_df_np_array = self._input_df_aqformat.to_numpy()

        self.extended_input_df = self._define_extended_input_df()

    
    def _define_extended_input_df(self):
        median_intensities_array = self._get_median_intensity_array()
        median_intensities_df = pd.DataFrame(median_intensities_array, index = self._input_df_aqformat.index, 
                                             columns = [f"median_rep{idx}" for idx in range(median_intensities_array.shape[1])])
        extended_input_df =  pd.concat([self._input_df_aqformat, median_intensities_df], axis="columns")
        extended_input_df = (2**extended_input_df).replace(np.nan, 0).reset_index()
        return extended_input_df
    
    def _get_median_intensity_array(self):
        per_condition_idxs = self._get_per_condition_idxs()
        max_num_replicates = max([len(x) for x in per_condition_idxs])

        median_intensities_array = np.full((self._input_df_np_array.shape[0], max_num_replicates), np.nan)

        for quant_id_idx in range(self._input_df_np_array.shape[0]):
            intensities = self._input_df_np_array[quant_id_idx,:]
            median_intensities = MedianIntensityCreator(intensities, 
                                                        per_condition_idxs, self._fraction_missing_values).median_intensities
            median_intensities_array[quant_id_idx, :len(median_intensities)] = median_intensities

        median_intensities_array = self._drop_allnan_columns(median_intensities_array)
        return median_intensities_array
    
    def _get_per_condition_idxs(self):
        samples_in_table = self._input_df_aqformat.columns.to_list()
        cond2samples = self._samplemap_df.groupby("condition")["sample"].apply(list).to_dict()
        idxs_per_condition = []
        for cond, samples in cond2samples.items():
            idxs = [samples_in_table.index(sample) for sample in samples]
            if len(idxs) >0:
                idxs_per_condition.append(idxs)
        
        return idxs_per_condition
    
    @staticmethod
    def _drop_allnan_columns(median_intensities_array):
        column_is_all_nan = np.all(np.isnan(median_intensities_array), axis = 0)
        return median_intensities_array[:, ~column_is_all_nan]
            
            

class MedianIntensityCreator():
    def __init__(self, intensities : np.array, per_condition_idxs : list[list], fraction_missing_values = 0.1):
        """This class takes the intensities for a given quant_id (length = samples*conditions).
        The median intensity is calculated by taking one intensity per condition (one replicate). As many median intensities 
        as possible are calculated obeying the fraction_missing_values parameter (dependent on the number of replicates with valid values).

        Args:
            intensities (np.array): intensities for all samples for a given quant_id. Corresponds to one row in the alphaquant input table.
            per_condition_idxs (list[list]): each list in the list contains the indices of the samples for a given condition.
            fraction_missing_values (float, optional): the maximum allowed fraction of missing values allowed when calcuating the median intensity. Defaults to 0.1.
        """
        self._intensities = intensities
        self._per_condition_idxs = per_condition_idxs
        self._fraction_missing_values = fraction_missing_values
        self._per_condition_intensities_sorted = self._define_per_condition_intensities_sorted()
        
        self.median_intensities = self._define_median_intensities()

    def _define_per_condition_intensities_sorted(self):
        per_condition_intensities = []
        for idxs in self._per_condition_idxs:
            condition_intensities = np.sort(self._intensities[idxs]) #nans are sorted to the end here
            per_condition_intensities.append(condition_intensities)
        return per_condition_intensities
    
    def _define_median_intensities(self):
        median_intensities = []

        max_num_replicates = max([len(x) for x in self._per_condition_idxs])
        for rep_idx in range(max_num_replicates):
            condition_intensities = self._get_intensities_for_replicate(rep_idx)
            fraction_missing = 1-len(condition_intensities)/max_num_replicates
            if fraction_missing <= self._fraction_missing_values:
                median_intensity = np.nanmedian(condition_intensities)
                median_intensities.append(median_intensity)
        return np.array(median_intensities)
    
    def _get_intensities_for_replicate(self, rep_idx):
         #rep_idx might be larger than the number of replicates for this particular condition as it goes until max_num_replicates
        condition_intensities = np.array([x[rep_idx] for x in self._per_condition_intensities_sorted 
                                              if self._check_if_replicate_is_valid(rep_idx, x)])
        return condition_intensities
    
    @staticmethod
    def _check_if_replicate_is_valid(rep_idx, condition_intensities):
        if rep_idx< len(condition_intensities):
            if ~np.isnan(condition_intensities[rep_idx]):
                return True
        return False


class ExtendedSampleMapCreator():
    def __init__(self, samplemap_df : pd.DataFrame, input_df_extended : pd.DataFrame):
        """This class extends the samplemap dataframe by adding the new median intensity samples

        Args:
            samplemap_df (pd.Dataframe): sample to condition mapping dataframe
            input_df_extended (pd.Dataframe): input dataframe with added median intensity samples
        """
        self._samplemap_df = samplemap_df
        self._input_df_extended = input_df_extended

        self.extended_samplemap_df = self._define_extended_samplemap_df()
    
    def _define_extended_samplemap_df(self):
        median_columns = [x for x in self._input_df_extended.columns if x.startswith("median_rep")]
        condition_name = ["median_reference" for x in range(len(median_columns))]
        median_name_df = pd.DataFrame({"sample": median_columns, "condition": condition_name})
        return pd.concat([self._samplemap_df, median_name_df], axis = "rows")