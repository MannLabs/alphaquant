import uuid
import pandas as pd
import numpy as np
import alphaquant.diffquant.background_distributions as aq_diff_bg
import pytest


def generate_random_input(num_pep, sample2cond_df):
    pepnames = list(map(lambda _idx: str(uuid.uuid4()), range(num_pep)))  # gives uuid strings for each peptide
    randarrays = 10 + 1.5 * np.random.randn(num_pep, sample2cond_df.shape[0])
    df_intens = pd.DataFrame(randarrays, columns=sample2cond_df["sample"].tolist())
    df_intens.insert(0, "peptides", pepnames)
    df_intens = df_intens.set_index("peptides")
    return df_intens

@pytest.fixture
def sample2cond_df():
    return pd.DataFrame({'sample': ['A1', 'A2', 'A3', 'B1', 'B2', 'B3'], 'condition': ['A', 'A', 'A', 'B', 'B', 'B']})

@pytest.fixture
def fixed_input(sample2cond_df):
    return generate_random_input(1000, sample2cond_df)

def test_condition_backgrounds(fixed_input):
    condbg = aq_diff_bg.ConditionBackgrounds(fixed_input, {})
    assert condbg.ion2background.keys() == condbg.ion2nonNanvals.keys()
