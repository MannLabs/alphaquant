import pytest
import pandas as pd
import numpy as np
import scipy.stats
import alphaquant.diffquant.background_distributions as aq_diff_bg
import alphaquant.diffquant.diff_analysis as aq_diff_diff


@pytest.fixture
def condbg():
    num_vals = 10000
    rand1 = np.random.normal(loc=0, size=num_vals)
    rand2 = np.random.normal(loc=0, size=num_vals)
    rand3 = np.random.normal(loc=0, size=num_vals)
    rand4 = np.random.normal(loc=0, size=num_vals)
    rand5 = np.random.normal(loc=0, size=num_vals)
    randarray = pd.DataFrame({1:rand1, 2:rand2, 3:rand3, 4:rand4, 5:rand5})
    return aq_diff_bg.ConditionBackgrounds(randarray, {})



def test_pval_assertion(condbg):
    noNanvals_from = np.array([1,2,1,2])
    noNanvals_to = np.array([2.1,1.1,1.9,0.9])
    pval = perform_peptide_difftest(condbg, noNanvals_from, 0, noNanvals_to, 1)
    assert pval > 0.95

def test_noreg_pepdist(condbg):
    num_samplings = 1000
    res = []
    for i in range(num_samplings):
        vals1 = 20 + np.random.randn(3)
        vals2 = 20 + np.random.randn(3)
        pval = perform_peptide_difftest(condbg, vals1, 0, vals2, 1)
        res.append(pval)
    
    unidist = np.random.uniform(size=1000)
    ks_stat, _ = scipy.stats.ks_2samp(res, unidist)
    assert ks_stat > 0.01

def perform_peptide_difftest(condbg, noNanvals_from, bg_idx1, noNanvals_to, bg_idx2):
    bgdist_from = condbg.backgrounds[bg_idx1]
    bgdist_to = condbg.backgrounds[bg_idx2]
    p2z = {}
    diffbg = aq_diff_bg.SubtractedBackgrounds(bgdist_from, bgdist_to, p2z)
    diffion = aq_diff_diff.DifferentialIon(noNanvals_from, noNanvals_to, diffbg, "", outlier_correction=False)
    return diffion.p_val