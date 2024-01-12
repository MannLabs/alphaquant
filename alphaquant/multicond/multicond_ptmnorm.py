import alphaquant.utils.utils as aqutils
import pandas as pd


def combine_results_tables(results_dir):
    condpair2results_df = MedianRefResultLoaderPTM(results_dir).condpair2results_df
    combiner = MedianConditionCombinerPTMNormed(condpair2results_df)
    combiner.combined_df.to_csv(f"{results_dir}/medianref_proteinnormed.tsv", sep = "\t")


class MedianRefResultLoaderPTM():
    def __init__(self, results_dir):
        self._results_dir = results_dir
        self.condpair2results_df = {}

        self._define_condpair2results_df()
    
    def _define_condpair2results_df(self):
        self.condpair2results_df = aqutils.read_all_results_files_in_results_folder(self._results_dir)
        self._remove_all_condpairs_w_no_medianref()
    
    def _remove_all_condpairs_w_no_medianref(self):
        keys_to_remove = [x for x in self.condpair2results_df.keys() if "median_reference" not in x]
        for key in keys_to_remove:
            del self.condpair2results_df[key]
    

class MedianConditionCombinerPTMNormed(): #this does a very simple combination of the conditionpairs, not taking into account proteoforms
    def __init__(self, condpair2results_df):
        self._condpair2results_df = condpair2results_df

        self.combined_df = None

        self._define_combined_df()

    def _define_combined_df(self):
        list_of_sub_dfs = []
        for condpair, results_df in self._condpair2results_df.items():
            sub_df = results_df.set_index("protein")["log2fc"].squeeze()
            sub_df.name = condpair
            list_of_sub_dfs.append(sub_df)
        self.combined_df = pd.concat(list_of_sub_dfs, axis=1)
    

