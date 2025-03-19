import alphaquant.multicond.median_condition_creation as aq_median_cond_creation
import numpy as np
import pandas as pd
import os
import tempfile


""" test MedianIntensityCreator """

def test_case_1():
    # No missing values, equal number of replicates per condition
    intensities = np.array([1, 1, 2, 2, 3, 3])
    per_condition_idxs = [[0, 1], [2, 3], [4, 5]]
    fraction_missing_values = 0.1
    
    mic = aq_median_cond_creation.MedianIntensityCreator(intensities, per_condition_idxs, fraction_missing_values)
    assert np.array_equal(mic.median_intensities, np.array([2, 2])), "Test Case 1 Failed"

def test_case_2():
    # No missing values, unequal number of replicates per condition v1
    intensities = np.array([1, 1, 2, 2,2, 3])
    per_condition_idxs = [[0, 1], [2, 3, 4], [5]]
    fraction_missing_values = 0.1
    
    mic = aq_median_cond_creation.MedianIntensityCreator(intensities, per_condition_idxs, fraction_missing_values)
    assert np.array_equal(mic.median_intensities, np.array([2])), "Test Case 2 Failed"


def test_case_3():
    # No missing values, unequal number of replicates per condition v2
    intensities = np.array([1, 1, 2, 2,2, 3, 3])
    per_condition_idxs = [[0, 1], [2, 3, 4], [5, 6]]
    fraction_missing_values = 0.1
    
    mic = aq_median_cond_creation.MedianIntensityCreator(intensities, per_condition_idxs, fraction_missing_values)
    assert np.array_equal(mic.median_intensities, np.array([2, 2])), "Test Case 3 Failed"

def test_case_4():
    # Missing values, unequal number of replicates per condition
    intensities = np.array([1, 1, 2, np.nan,np.nan, 3, 3])
    per_condition_idxs = [[0, 1], [2, 3, 4], [5, 6]]
    fraction_missing_values = 0.1
    
    mic = aq_median_cond_creation.MedianIntensityCreator(intensities, per_condition_idxs, fraction_missing_values)
    assert np.array_equal(mic.median_intensities, np.array([2])), "Test Case 4 Failed"

def test_case_5():
    # Missing values v3, unequal number of replicates per condition
    intensities = np.array([np.nan, 1, 2, np.nan,np.nan, np.nan, 3])
    per_condition_idxs = [[0, 1], [2, 3, 4], [5, 6]]
    fraction_missing_values = 0.1
    
    mic = aq_median_cond_creation.MedianIntensityCreator(intensities, per_condition_idxs, fraction_missing_values)
    assert np.array_equal(mic.median_intensities, np.array([2])), "Test Case 5 Failed"


def test_case_6():
    # Edge case: fraction_missing_values = 1
    intensities = np.array([1, 1, 2, np.nan,np.nan, 3, 3])
    per_condition_idxs = [[0, 1], [2, 3, 4], [5, 6]]
    fraction_missing_values = 1
    
    mic = aq_median_cond_creation.MedianIntensityCreator(intensities, per_condition_idxs, fraction_missing_values)
    assert np.allclose(mic.median_intensities, np.array([2, 2, np.nan]), equal_nan = True), "Test Case 6 Failed"

def test_case_7():
    # Edge case: fraction_missing_values = 0
    intensities = np.array([1, 1, 2, np.nan,np.nan, 3, 3])
    per_condition_idxs = [[0, 1], [2, 3, 4], [5, 6]]
    fraction_missing_values = 0
    
    mic = aq_median_cond_creation.MedianIntensityCreator(intensities, per_condition_idxs, fraction_missing_values)
    assert np.array_equal(mic.median_intensities, np.array([2]), equal_nan = True), "Test Case 7 Failed"

def test_case_8():
    # Edge case: all nan values
    intensities = np.array([np.nan, np.nan, np.nan, np.nan,np.nan, np.nan, np.nan])
    per_condition_idxs = [[0, 1], [2, 3, 4], [5, 6]]
    fraction_missing_values = 0
    
    mic = aq_median_cond_creation.MedianIntensityCreator(intensities, per_condition_idxs, fraction_missing_values)
    assert np.array_equal(mic.median_intensities, np.array([])), "Test Case 8 Failed"





"""test ExtendedInputDfCreator"""


def test_median_condition_creator():

    test_intensities = [
    [1, 1, 2, 0, 2, 3, 3],
    [1, 1, 2, 2, 2, 0, 3],
    [1, 1, 2, 2, 2, 3, 3],
    [1, 1, 2, 0, 0, 3, 3],
    [0, 1, 2, 0, 0, 0, 3],
    ]

    test_intensities_2 = [np.array(row) *2  for row in test_intensities]
    test_intensities = test_intensities + test_intensities_2

    input_df_test = pd.DataFrame(test_intensities)



    sample_names = [f"sample_{x}" for x in range(len(input_df_test.columns))]
    condition_names = ["condition_1", "condition_1", "condition_2", "condition_2", "condition_2", "condition_3", "condition_3"]
    samplemap_df = pd.DataFrame({"sample": sample_names, "condition": condition_names})

    input_df_test.columns = sample_names
    input_df_test["quant_id"] = [f"quant_{x}" for x in range(len(input_df_test.index))]
    input_df_test["protein"] = ["protein_1", "protein_1", "protein_1", "protein_1", "protein_1", "protein_1", "protein_2", "protein_3", "protein_3", "protein_3"]
    mic = aq_median_cond_creation.ExtendedInputDfCreator(input_df_test, samplemap_df, fraction_missing_values = 0.1)

    assert np.array_equal(mic.extended_input_df.loc[:,"median_rep0"], np.array([2, 2, 2, 2, 2, 4, 4, 4, 4, 4])), "Test Case Failed"

    print("Test Passed")


"""test MedianConditionManager"""



def generate_samplemap_df(replicate_numbers):
        
    num_conditions = len(replicate_numbers)
    samplemap_data = {
        'sample': [],
        'condition': []
    }
    
    for condition_idx in range(num_conditions):
        condition_name = f'condition{condition_idx + 1}'
        for replicate_idx in range(replicate_numbers[condition_idx]):
            samplemap_data['sample'].append(f'sample{condition_idx + 1}_{replicate_idx + 1}')
            samplemap_data['condition'].append(condition_name)
            
    return pd.DataFrame(samplemap_data)

def generate_normalized_input_df(samplemap_df, num_quant_ids =10, nan_fraction = 0.1, add_additional_samples = False):
    np.random.seed(0)  # For reproducibility
    quant_ids = [f"quant_{i}" for i in range(num_quant_ids)]
    proteins = [f"protein_{i}" for i in range(num_quant_ids)]
    reformatted_input_data = {'quant_id': quant_ids, 'protein': proteins}

    sample_list = get_sample_list(samplemap_df, add_additional_samples)
    
    for sample in sample_list:
        data = np.random.rand(num_quant_ids)
        if nan_fraction > 0:


            nan_indices = np.random.choice(num_quant_ids, int(np.floor(nan_fraction * num_quant_ids)), replace=False)
            data[nan_indices] = np.nan
        reformatted_input_data[sample] = data
        
    return pd.DataFrame(reformatted_input_data)

def get_sample_list(samplemap_df, add_additional_samples):
    sample_list = samplemap_df['sample'].to_list()
    if add_additional_samples:
        sample_list.extend(['additional_sample1', 'additional_sample2'])
    return sample_list



def test_add_and_save_median_condition():
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create temporary file paths
        input_file_path = os.path.join(temp_dir, 'input.tsv')
        samplemap_file_path = os.path.join(temp_dir, 'samplemap.tsv')

        # Generate test data and save to temporary files
        replicate_numbers = [3, 2, 3]
        samplemap_df = generate_samplemap_df(replicate_numbers)
        samplemap_df.to_csv(samplemap_file_path, sep="\t", index=False)

        input_df = generate_normalized_input_df(samplemap_df, num_quant_ids=100, nan_fraction=0.3)
        input_df.to_csv(input_file_path, sep="\t", index=False)

        # Call function under test
        median_condition_manager = aq_median_cond_creation.MedianConditionManager(input_file_path, samplemap_file_path)
        
        input_file_path_adapted = median_condition_manager.input_filename_adapted
        samplemap_file_path_adapted = median_condition_manager.samplemap_filename_adapted

        # Check that the files were created
        assert os.path.exists(input_file_path_adapted), "Input file was not created."
        assert os.path.exists(samplemap_file_path_adapted), "Samplemap file was not created."

        # Read the files back in and validate contents
        extended_input_df_saved = pd.read_csv(input_file_path_adapted, sep="\t")
        extended_samplemap_df_saved = pd.read_csv(samplemap_file_path_adapted, sep="\t")

        # Validate that the DataFrames have the added median conditions
        median_condition_creator = aq_median_cond_creation.ExtendedInputDfCreator(input_df, samplemap_df)
        pd.testing.assert_frame_equal(extended_input_df_saved, median_condition_creator.extended_input_df)     
        pd.testing.assert_frame_equal(extended_samplemap_df_saved.reset_index(drop= True), 
                                      median_condition_manager.samplemap_df_extended.reset_index(drop = True))

        print("Test for add_and_save_median_condition passed.")


