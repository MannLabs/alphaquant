import numpy as np 
import pandas as pd
import alphaquant.diffquant.condpair_analysis as aq_diffquant_condpair


def test_get_per_condition_dataframes():
        data = {'A': [1, np.nan, np.nan, 4, 5, 123, 3, 4],
        'B': [5, np.nan, np.nan, 8, 9, 123, 4, 5],
        'C': [10, 11, 5, np.nan, 14, 234, 5, 6],
        'D': [15, 16, 6, np.nan, np.nan, 3, 55, 6]}
        unnormed_df = pd.DataFrame(data)

        # Test 1
        samples_c1 = ['A', 'B']
        samples_c2 = ['C', 'D']
        minrep_both = 2
        minrep_either = None
        minrep_c1 = 2 
        minrep_c2 = 0 #minrep c1 and c2 should override minrep both which is set per default
        df_c1, df_c2 = aq_diffquant_condpair.get_per_condition_dataframes(samples_c1, samples_c2, unnormed_df, minrep_both, minrep_either, minrep_c1, minrep_c2)

        assert not df_c1.empty and not df_c2.empty, "DataFrames should not be empty"
        assert len(df_c1.index) >= 5 and len(df_c2.index) >= 5, "DataFrames should have enough data"
        assert get_minrep_of_df(df_c1) >= 2, "DataFrame should have minrep of 2"
        assert get_minrep_of_df(df_c2) == 0, "DataFrame should have minrep of 0"

        # Test 2
        samples_c1 = ['A', 'B']
        samples_c2 = ['C', 'D']
        minrep_both = 2
        minrep_either = None
        minrep_c1 = None
        minrep_c2 = None

        df_c1, df_c2 = aq_diffquant_condpair.get_per_condition_dataframes(samples_c1, samples_c2, unnormed_df, minrep_both, minrep_either, minrep_c1, minrep_c2)

        assert not df_c1.empty and not df_c2.empty, "DataFrames should not be empty"
        assert len(df_c1.index) >= 5 and len(df_c2.index) >= 5, "DataFrames should have enough data"
        assert get_minrep_of_df(df_c1) >= 2 and get_minrep_of_df(df_c2) >= 2, "DataFrames should have minrep of 2"

        # Test 3
        samples_c1 = ['A', 'B']
        samples_c2 = ['C', 'D']
        minrep_both = None
        minrep_either = 2
        minrep_c1 = None
        minrep_c2 = None

        df_c1, df_c2 = aq_diffquant_condpair.get_per_condition_dataframes(samples_c1, samples_c2, unnormed_df, minrep_both, minrep_either, minrep_c1, minrep_c2)

        assert not df_c1.empty and not df_c2.empty, "DataFrames should not be empty"
        assert len(df_c1.index) >= 5 and len(df_c2.index) >= 5, "DataFrames should have enough data"
        assert get_minrep_of_df(df_c1) ==0  and get_minrep_of_df(df_c2)  == 0, "DataFrames should have minrep of 0"


        # Test 4
        samples_c1 = ['A', 'B']
        samples_c2 = ['C', 'D']
        minrep_both = 2
        minrep_either = 2 #minrep either should override minrep both which is set per default
        minrep_c1 = None
        minrep_c2 = None

        df_c1, df_c2 = aq_diffquant_condpair.get_per_condition_dataframes(samples_c1, samples_c2, unnormed_df, minrep_both, minrep_either, minrep_c1, minrep_c2)

        assert not df_c1.empty and not df_c2.empty, "DataFrames should not be empty"
        assert len(df_c1.index) >= 5 and len(df_c2.index) >= 5, "DataFrames should have enough data"
        assert get_minrep_of_df(df_c1) ==0  and get_minrep_of_df(df_c2)  == 0, "DataFrames should have minrep of 0"



        print("All tests passed")



def get_minrep_of_df(df):
    num_not_nan_per_row = df.notna().sum(axis=1)
    return num_not_nan_per_row.min()




