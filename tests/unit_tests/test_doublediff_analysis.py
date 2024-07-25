import pandas as pd
import numpy as np
import scipy.stats as stats
import alphaquant.diffquant.doublediff_analysis as aq_diff_doublediff
import alphaquant.diffquant.background_distributions as aq_diff_background
import random

P2Z = {}
ION2DIFFDIST = {}
DEEDPAIR2DOUBLEDIFFDIST = {}



def test_that_pvalues_are_uniformly_distributed_multiple_scenarios():
    # Test case 1: Most difficult case with different distributions, protein-specific shift, and NaN values
    check_pvalues_uniform_distribution(
        num_samples_c1=3, num_samples_c2=5, num_proteins=200,
        c1_distribution_width=50, c2_distribution_width=10,
        add_protein_specific_shift=True,
        c1_nan_fraction=0.2, c2_nan_fraction=0.4,
        seed=42
    )

    # Test case 2: Equal samples, no shift
    check_pvalues_uniform_distribution(
        num_samples_c1=5, num_samples_c2=5, num_proteins=100,
        c1_distribution_width=30, c2_distribution_width=30,
        add_protein_specific_shift=False,
        c1_nan_fraction=0.1, c2_nan_fraction=0.1,
        seed=42
    )

    # Test case 3: Unequal samples, with shift
    check_pvalues_uniform_distribution(
        num_samples_c1=4, num_samples_c2=6, num_proteins=75,
        c1_distribution_width=40, c2_distribution_width=20,
        add_protein_specific_shift=True,
        c1_nan_fraction=0.3, c2_nan_fraction=0.3,
        seed=42
    )

    # Test case 4: Reverse unequal samples, no shift
    check_pvalues_uniform_distribution(
        num_samples_c1=6, num_samples_c2=4, num_proteins=60,
        c1_distribution_width=25, c2_distribution_width=25,
        add_protein_specific_shift=False,
        c1_nan_fraction=0.15, c2_nan_fraction=0.25,
        seed=42
    )


def check_pvalues_uniform_distribution(num_samples_c1, num_samples_c2, num_proteins, 
                                       c1_distribution_width, c2_distribution_width, 
                                       add_protein_specific_shift, c1_nan_fraction, 
                                       c2_nan_fraction, seed):
    """
    This function checks whether the p-values from the double differential test are uniformly distributed as expected.
    We test the "most difficult" case, with different distributions for the two conditions, a protein-specific shift,
    Nan values in the data, and a different number of samples in the two conditions. 
    """

    simulated_df = _simulate_alphaquant_input_table(
        num_samples_c1=num_samples_c1, 
        num_samples_c2=num_samples_c2, 
        num_proteins=num_proteins, 
        c1_distribution_width=c1_distribution_width, 
        c2_distribution_width=c2_distribution_width, 
        add_protein_specific_shift=add_protein_specific_shift, 
        c1_nan_fraction=c1_nan_fraction, 
        c2_nan_fraction=c2_nan_fraction, 
        seed=seed
    )
    samplemap_df = _create_samplemap_from_df_simulated(simulated_df)
    pvals = _compare_sets_of_peptides(simulated_df, samplemap_df)

    np.random.seed(42)
    unidist = np.random.uniform(size=1000)

    pval_ks = stats.ks_2samp(pvals, unidist)[1]
    assert pval_ks > 0.05, f"P-value {pval_ks} is not greater than 0.05"
    print(f"Test passed. KS test shows no rejection of null hypothesis as expected.")





def _simulate_alphaquant_input_table(num_samples_c1, num_samples_c2, num_proteins, 
                                    peptides_per_protein_range=(1, 7),
                                    c1_distribution_width=20,
                                    c2_distribution_width=20,
                                    add_protein_specific_shift=False,
                                    c1_nan_fraction=0,
                                    c2_nan_fraction=0,
                                    seed=42):
    """
    This function simulates an input table for the alphaquant package and returns it as a pandas DataFrame. The table will 
    simulate peptides for a given protein all of which are sampled from the same distribution. This way, performing a
    double differential test by comparing two sets of peptides from the same protein should yield a uniform distribution
    of p-values.
    
    Args:
    num_samples_c1 (int): Number of samples in condition 1
    num_samples_c2 (int): Number of samples in condition 2
    num_proteins (int): Number of proteins to simulate
    peptides_per_protein_range (tuple, optional): Range of peptides per protein. Defaults to (1, 7).
    c1_distribution_width (float): Width (standard deviation) of the distribution for condition 1. Defaults to 20.
    c2_distribution_width (float): Width (standard deviation) of the distribution for condition 2. Defaults to 20.
    add_protein_specific_shift (bool): If True, adds a protein-specific shift to condition 2. Defaults to False.
    c1_nan_fraction (float): Fraction of values to be replaced with NaN in condition 1. Defaults to 0.
    c2_nan_fraction (float): Fraction of values to be replaced with NaN in condition 2. Defaults to 0.
    
    Returns:
    pandas.DataFrame: Simulated input table for alphaquant
    """
    
    np.random.seed(seed)
    # Initialize empty lists to store data
    proteins = []
    quant_ids = []
    sample_data = []
    
    # Generate data for each protein
    for i in range(1, num_proteins + 1):
        protein = f"protein_{i}"
        num_peptides = np.random.randint(peptides_per_protein_range[0], peptides_per_protein_range[1] + 1)
        
        # Generate protein-specific shift if required
        protein_shift = np.random.uniform(50, 150) if add_protein_specific_shift else 0
        
        for j in range(1, num_peptides + 1):
            proteins.append(protein)
            quant_ids.append(f"{protein}_peptide_{j}")
            
            # Generate sample data for condition 1
            c1_values = np.abs(np.random.normal(loc=100, scale=c1_distribution_width, size=num_samples_c1))
            
            # Generate sample data for condition 2
            c2_values = np.abs(np.random.normal(loc=100, scale=c2_distribution_width, size=num_samples_c2))
            
            # Add protein-specific shift to condition 2 if required
            if add_protein_specific_shift:
                c2_values += protein_shift
            
            # Introduce NaN values while ensuring at least one valid value per condition
            c1_values = _introduce_nans(c1_values, c1_nan_fraction)
            c2_values = _introduce_nans(c2_values, c2_nan_fraction)
            
            # Combine the values
            sample_values = np.concatenate([c1_values, c2_values])
            sample_data.append(sample_values)
    
    # Create DataFrame
    simulated_df = pd.DataFrame({
        'protein': proteins,
        'quant_id': quant_ids
    })
    
    # Add sample columns
    for i in range(1, num_samples_c1 + 1):
        simulated_df[f'c1_sample_{i}'] = [row[i-1] for row in sample_data]
    
    for i in range(1, num_samples_c2 + 1):
        simulated_df[f'c2_sample_{i}'] = [row[i-1+num_samples_c1] for row in sample_data]
    
    simulated_df = simulated_df.set_index('quant_id')
    return simulated_df

def _introduce_nans(values, nan_fraction):
        """Helper function to introduce NaNs while ensuring at least one valid value"""
        n = len(values)
        nan_count = int(n * nan_fraction)
        if nan_count >= n:
            nan_count = n - 1  # Ensure at least one value remains
        nan_indices = np.random.choice(n, nan_count, replace=False)
        values[nan_indices] = np.nan
        return values


def _create_samplemap_from_df_simulated(simulated_df):
    """
    This function reads in the columns of the simulated input table and creates a dataframe 
    with the following columns:
    - sample (sample name)
    - condition (condition name)
    
    Args:
    simulated_df (pandas.DataFrame): The simulated input table DataFrame
    
    Returns:
    pandas.DataFrame: A DataFrame containing sample names and their corresponding conditions
    """
    
    # Get all column names
    all_columns = simulated_df.columns
    
    # Filter for sample columns (those starting with 'c1_' or 'c2_')
    sample_columns = [col for col in all_columns if col.startswith(('c1_', 'c2_'))]
    
    # Create lists to store sample names and conditions
    samples = []
    conditions = []
    
    # Populate the lists
    for column in sample_columns:
        samples.append(column)
        condition = 'condition_1' if column.startswith('c1_') else 'condition_2'
        conditions.append(condition)
    
    # Create the DataFrame
    samplemap_df = pd.DataFrame({
        'sample': samples,
        'condition': conditions
    })
    
    return samplemap_df


    
def _compare_sets_of_peptides(simulated_df, samplemap_df):
    """
    This function iterates through every protein in the simulated_df and randomly subsets the peptides into two sets. 
    It then calculates the doublediff score for the two sets of peptides.

    Args:
    simulated_df (pandas.DataFrame): The simulated input DataFrame
    samplemap_df (pandas.DataFrame): The sample map DataFrame

    Returns:
    dict: A dictionary with protein names as keys and their doublediff scores as values
    """
    normed_c1, normed_c2 = _initialize_condition_backgrounds(simulated_df, samplemap_df)
    
    pvals = []
    
    for protein in simulated_df['protein'].unique():
        protein_peptides = simulated_df[simulated_df['protein'] == protein].index.tolist()
        
        # Randomly split peptides into two sets
        n_peptides = len(protein_peptides)
        if n_peptides < 2:
            continue  # Skip proteins with less than 2 peptides
        
        split_point = random.randint(1, n_peptides - 1)
        ions1 = protein_peptides[:split_point]
        ions2 = protein_peptides[split_point:]
        
        # Calculate doublediff score

        fcfc, p = aq_diff_doublediff.calc_doublediff_score(ions1=ions1, ions2=ions2, normed_c1=normed_c1, normed_c2=normed_c2, ion2diffDist=ION2DIFFDIST, p2z=P2Z, deedpair2doublediffdist=DEEDPAIR2DOUBLEDIFFDIST)
        pvals.append(p)
    
    return pvals


def _initialize_condition_backgrounds(simulated_df, samplemap_df):
    df_c1, df_c2 = _get_c1_c2_dfs(simulated_df, samplemap_df, ["condition_1", "condition_2"])
    normed_c1 = aq_diff_background.ConditionBackgrounds(df_c1, P2Z)
    normed_c2 = aq_diff_background.ConditionBackgrounds(df_c2, P2Z)

    return normed_c1, normed_c2

def _get_c1_c2_dfs(unnormed_df, labelmap_df, condpair, minrep = 2):
    c1_samples = labelmap_df[labelmap_df["condition"]== condpair[0]]
    c2_samples = labelmap_df[labelmap_df["condition"]== condpair[1]]
    df_c1 = unnormed_df.loc[:, c1_samples["sample"]].dropna(thresh=minrep, axis=0)
    df_c2 = unnormed_df.loc[:, c2_samples["sample"]].dropna(thresh=minrep, axis=0)

    return df_c1, df_c2



