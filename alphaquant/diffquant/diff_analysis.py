
from statistics import NormalDist
import numpy as np
import math
import statistics
from scipy.stats import ttest_ind
from scipy.stats import t as student_t
import alphaquant.diffquant.diffutils as aqdiffutils

class DifferentialIon():
    """Computes differential statistics for an ion using empirical background distributions.

    This is the default statistical test in AlphaQuant. It uses intensity-dependent empirical
    background distributions to compute p-values and z-scores for each ion. The method accounts
    for technical variation by comparing observed fold changes against distributions derived from
    similarly abundant ions in the dataset.

    Args:
        noNanvals_from: Array of log2 intensities from condition 1 (no NaN values)
        noNanvals_to: Array of log2 intensities from condition 2 (no NaN values)
        diffDist: BackGroundDistribution object containing empirical error distribution
        name: Ion identifier (e.g., peptide sequence + charge state)
        outlier_correction: If True, inflates variance estimates when replicates show
                          unusually high variability (default: True)

    Attributes:
        name: Ion identifier
        p_val: Two-sided p-value testing null hypothesis of no change
        fc: Log2 fold change (condition1 - condition2)
        z_val: Z-score (signed, magnitude indicates significance)
        usable: Boolean indicating whether statistics could be computed
    """

    def __init__(self,noNanvals_from, noNanvals_to, diffDist, name, outlier_correction = True):


        self.name = name
        self.p_val = None
        self.fc = None
        self.z_val = None
        self.usable = False

        self._calc_diffreg_peptide(noNanvals_from, noNanvals_to, diffDist, outlier_correction)


    def _calc_diffreg_peptide(self, noNanvals_from, noNanvals_to, diffDist, outlier_correction):
        """Calculates differential expression statistics using empirical background distributions.

        This method computes all pairwise fold changes between replicates of the two conditions,
        converts them to z-scores using the empirical background distribution, and aggregates
        them into a single p-value. The variance is scaled based on the number of replicates
        and optionally adjusted for outliers.

        Args:
            noNanvals_from: Log2 intensities from condition 1
            noNanvals_to: Log2 intensities from condition 2
            diffDist: Empirical background distribution for variance estimation
            outlier_correction: Whether to apply robust variance inflation
        """

        nrep_from = len(noNanvals_from)
        nrep_to = len(noNanvals_to)

        if ((nrep_from==0) or (nrep_to ==0)):
            return
        var_from = diffDist.var_from
        var_to = diffDist.var_to

        perEvidenceVariance = diffDist.var + (nrep_to-1) * var_from + (nrep_from-1) * var_to
        totalVariance = perEvidenceVariance*nrep_to * nrep_from
        outlier_scaling_factor = 1.0
        if outlier_correction:
            outlier_scaling_factor = calc_outlier_scaling_factor(noNanvals_from, noNanvals_to, diffDist)

        fc_sum =0
        z_sum=0
        unscaled_zs = []
        for from_intens in noNanvals_from:
            for to_intens in noNanvals_to:
                fc = from_intens - to_intens
                fc_sum+=fc
                z_unscaled = diffDist.calc_zscore_from_fc(fc)
                unscaled_zs.append(z_unscaled)
                z_sum += z_unscaled/outlier_scaling_factor

        scaled_SD =  math.sqrt(totalVariance/diffDist.var)*outlier_scaling_factor

        self.fc = fc_sum/(nrep_from * nrep_to)
        self.p_val = 2.0 * (1.0 -  NormalDist(mu=0, sigma= scaled_SD).cdf(abs(z_sum)))
        self.z_val = z_sum/scaled_SD
        self.usable = True




  #self.var_from = from_dist.var
   #     self.var_to

class DifferentialIonTTest():
    """Example implementation of differential testing using Welch's t-test.

    This is an alternative statistical test provided as example code to demonstrate how
    developers can implement their own methods. It uses Welch's two-sample t-test with
    robust variance estimation (similar to MS-EmpiRe). This implementation has not been
    extensively benchmarked and is included primarily for educational purposes.

    To use this method instead of DifferentialIon, set ion_test_method='ttest' when
    calling run_pipeline().

    Args:
        noNanvals_from: Array of log2 intensities from condition 1 (no NaN values)
        noNanvals_to: Array of log2 intensities from condition 2 (no NaN values)
        name: Ion identifier (e.g., peptide sequence + charge state)
        p2z: Dictionary cache for p-value to z-value conversions
        outlier_correction: If True, inflates standard error using robust estimators
                          (MAD/IQR) when sample variance is unusually low (default: True)

    Attributes:
        name: Ion identifier
        p_val: Two-sided p-value from Welch's t-test
        fc: Log2 fold change (condition1 - condition2)
        z_val: Z-score derived from p-value
        usable: Boolean indicating whether statistics could be computed
    """

    def __init__(self, noNanvals_from, noNanvals_to, name, p2z = None, outlier_correction: bool = True):


        self.name = name
        self.p_val = None
        self.fc = None
        self.z_val = None
        self.usable = False

        self._calc_ttest_peptide(noNanvals_from, noNanvals_to, p2z, outlier_correction)


    def _calc_ttest_peptide(self, noNanvals_from, noNanvals_to, p2z, outlier_correction):
        """Calculates differential expression statistics using Welch's t-test.

        Computes the t-statistic and p-value using scipy's Welch's t-test, with optional
        robust variance inflation. The p-value is then converted to a z-score for
        compatibility with the tree aggregation framework.

        Args:
            noNanvals_from: Log2 intensities from condition 1
            noNanvals_to: Log2 intensities from condition 2
            p2z: Cache dictionary for p-value to z-value conversions
            outlier_correction: Whether to apply robust standard error inflation
        """

        nrep_from = len(noNanvals_from)
        nrep_to = len(noNanvals_to)

        if ((nrep_from==0) or (nrep_to ==0)):
            return

        mean_from = float(np.mean(noNanvals_from))
        mean_to = float(np.mean(noNanvals_to))
        self.fc = mean_from - mean_to

        # Compute Welch's t-statistic and degrees of freedom
        n1 = nrep_from
        n2 = nrep_to
        try:
            res = ttest_ind(noNanvals_from, noNanvals_to, equal_var=False, nan_policy='omit')
            t_stat = float(res.statistic) if res.statistic is not None else 0.0
            p_val = float(res.pvalue) if res.pvalue is not None else 1.0
        except Exception:
            t_stat = 0.0
            p_val = 1.0

        s1 = float(np.std(noNanvals_from, ddof=1)) if n1 > 1 else 0.0
        s2 = float(np.std(noNanvals_to, ddof=1)) if n2 > 1 else 0.0
        se_standard_sq = (s1*s1)/max(1, n1) + (s2*s2)/max(1, n2)
        se_standard = math.sqrt(se_standard_sq) if se_standard_sq > 0 else 0.0

        # Welch-Satterthwaite degrees of freedom
        num_df = se_standard_sq * se_standard_sq
        den_df = 0.0
        if n1 > 1:
            den_df += ((s1*s1)/n1) * ((s1*s1)/n1) / (n1 - 1)
        if n2 > 1:
            den_df += ((s2*s2)/n2) * ((s2*s2)/n2) / (n2 - 1)
        df = (num_df/den_df) if den_df > 0 else max(n1 + n2 - 2, 1)

        # Robust SE inflation analogous to diffdist outlier scaling
        if outlier_correction and se_standard > 0 and n1 > 1 and n2 > 1:
            se_robust = _calc_robust_se_ttest(noNanvals_from, noNanvals_to)
            if se_robust > 0:
                scaling_factor = max(1.0, min(5.0, se_robust / se_standard))
                t_adj = t_stat / scaling_factor
                p_val = 2.0 * float(student_t.sf(abs(t_adj), df))

        abs_z = aqdiffutils.two_sided_p_to_abs_z(p_val, p2z if p2z is not None else {})
        sign = 1.0 if self.fc >= 0 else -1.0
        self.z_val = sign * abs_z
        self.p_val = p_val
        self.usable = True


def calc_outlier_scaling_factor(noNanvals_from, noNanvals_to, diffDist):
    """Computes a variance inflation factor to account for outlier replicates.

    Compares the between-replicate variance to the expected technical variance from the
    empirical background distribution. If replicates are more variable than expected
    (e.g., due to biological variability or technical outliers), the variance estimate
    is inflated accordingly. This makes the test more conservative when data quality is poor.

    Args:
        noNanvals_from: Log2 intensities from condition 1
        noNanvals_to: Log2 intensities from condition 2
        diffDist: Background distribution providing expected technical variance

    Returns:
        float: Scaling factor (>=1.0) to multiply variance by. Returns 1.0 if no
               correction needed, higher values when outliers are detected.
    """
    sd_from = math.sqrt(diffDist.var_from)
    sd_to = math.sqrt(diffDist.var_to)
    median_from = statistics.median(noNanvals_from)
    median_to = statistics.median(noNanvals_to)

    between_rep_SD_from = math.sqrt(sum(np.square(noNanvals_from-median_from))/len(noNanvals_from)) if len(noNanvals_from)>1 else sd_from
    between_rep_SD_to = math.sqrt(sum(np.square(noNanvals_to-median_to))/len(noNanvals_to)) if len(noNanvals_to)>1 else sd_to

    highest_SD_from = max(between_rep_SD_from, sd_from)
    highest_SD_to = max(between_rep_SD_to, sd_to)
    highest_SD_combined = math.sqrt(highest_SD_from**2 + highest_SD_to**2)

    scaling_factor = max(1.0, highest_SD_combined/diffDist.SD)
    return scaling_factor

def _robust_sd(x):
    """Computes a robust estimate of standard deviation using MAD or IQR.

    This function provides a robust alternative to standard deviation that is less
    sensitive to outliers. It tries three methods in order of robustness:
    1. MAD (Median Absolute Deviation) scaled to match SD for normal distributions
    2. IQR (Interquartile Range) scaled to match SD for normal distributions
    3. Regular sample standard deviation (fallback)

    Args:
        x: Array of values

    Returns:
        float: Robust estimate of standard deviation, or 0.0 if n < 2
    """
    x = np.asarray(x)
    n = x.size
    if n == 0:
        return 0.0
    if n == 1:
        return 0.0
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    if mad > 0:
        return 1.4826 * mad
    # Fallback to IQR if MAD is zero
    q75, q25 = np.percentile(x, [75, 25])
    iqr = float(q75 - q25)
    if iqr > 0:
        return iqr / 1.349
    # Final fallback to sample SD
    return float(np.std(x, ddof=1))

def _calc_robust_se_ttest(noNanvals_from, noNanvals_to):
    """Computes robust standard error for t-test using inflated variance estimates.

    This function calculates the standard error for Welch's t-test, but uses the maximum
    of the regular standard deviation and a robust estimate (MAD/IQR-based). This provides
    protection against underestimation of variance when sample SD is unusually low.

    Args:
        noNanvals_from: Log2 intensities from condition 1
        noNanvals_to: Log2 intensities from condition 2

    Returns:
        float: Robust standard error estimate, or 0.0 if insufficient replicates
    """
    n1 = len(noNanvals_from)
    n2 = len(noNanvals_to)
    if n1 < 2 or n2 < 2:
        return 0.0
    s1 = float(np.std(noNanvals_from, ddof=1))
    s2 = float(np.std(noNanvals_to, ddof=1))
    s1_rob = _robust_sd(noNanvals_from)
    s2_rob = _robust_sd(noNanvals_to)
    s1_infl = max(s1, s1_rob)
    s2_infl = max(s2, s2_rob)
    se_sq = (s1_infl*s1_infl)/n1 + (s2_infl*s2_infl)/n2
    return math.sqrt(se_sq) if se_sq > 0 else 0.0
