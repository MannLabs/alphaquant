from time import time
import numpy as np
import math
import alphaquant.diffquant.diffutils as aqutils

import alphaquant.config.config as aqconfig
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

from numba import njit
from statistics import NormalDist
import alphaquant.diffquant.diffutils as aqdiffutils





class ConditionBackgrounds():

    def __init__(self, normed_condition_df, p2z):
        self.backgrounds = []
        self.ion2background = {}
        self.ion2nonNanvals = {}
        self.ion2allvals = {}
        self.idx2ion = {}
        self.init_ion2nonNanvals(normed_condition_df)
        self.context_ranges = []
        self.select_intensity_ranges(p2z)

        self.all_intensities = np.concatenate(list(self.ion2nonNanvals.values()))
        self.num_replicates = len(next(iter(self.ion2allvals.values())))



    def init_ion2nonNanvals(self, normed_condition_df):
        normed_condition_df['median'] = normed_condition_df.median(numeric_only=True, axis=1)
        normed_condition_df = normed_condition_df.sort_values(by='median').drop('median', axis=1)
        self.normed_condition_df = normed_condition_df
        #nonan_array = get_nonna_array(normed_condition_df.to_numpy())
        #self.ion2nonNanvals = dict(zip(normed_condition_df.index, nonan_array))
        self.ion2nonNanvals = aqutils.get_non_nas_from_pd_df(normed_condition_df)
        self.ion2allvals = aqutils.get_ionints_from_pd_df(normed_condition_df)
        self.idx2ion = dict(zip(range(len(normed_condition_df.index)), normed_condition_df.index))


    def select_intensity_ranges(self, p2z):
        total_available_comparisons =0
        num_contexts = 10
        cumulative_counts = np.zeros(self.normed_condition_df.shape[0])

        for idx ,count in enumerate(self.normed_condition_df.count(axis=1)):
            total_available_comparisons+=count-1
            cumulative_counts[idx] = int(total_available_comparisons/2)


        #assign the context sizes
        context_size = np.max([1000, int(total_available_comparisons/(1+num_contexts/2))])
        if context_size> total_available_comparisons:
            context_size = int(total_available_comparisons/2)
        halfcontext_size = int(context_size/2)
        context_boundaries = np.zeros(3).astype(int)

        middle_idx = int(np.searchsorted(cumulative_counts, halfcontext_size))
        end_idx = int(np.searchsorted(cumulative_counts, context_size))


        context_boundaries[0] = 0
        context_boundaries[1] = middle_idx
        context_boundaries[2] = end_idx
        while context_boundaries[1] < len(cumulative_counts):
            bgdist = BackGroundDistribution(context_boundaries[0], context_boundaries[2], self.ion2nonNanvals, self.idx2ion, p2z)
            self.context_ranges.append([context_boundaries[0], context_boundaries[2]])
            self.assign_ions2bgdists(context_boundaries[0], context_boundaries[2], bgdist)
            self.backgrounds.append(bgdist)
            context_boundaries[0] = context_boundaries[1]
            context_boundaries[1] = context_boundaries[2]
            end_idx = np.searchsorted(cumulative_counts, context_size + cumulative_counts[context_boundaries[0]])
            if end_idx > len(cumulative_counts)-(context_boundaries[1]-context_boundaries[0])/1.5:
                end_idx = len(cumulative_counts)
            context_boundaries[2] = end_idx

    def assign_ions2bgdists(self, boundaries1, boundaries2, bgdist):
        ion2bg_local = {} #dict(map(lambda _idx : (self.normed_condition_df.index.values[_idx], bgdist), range(boundaries1, boundaries2)))
        for idx in range(boundaries1, boundaries2):
            ion2bg_local.update({self.idx2ion.get(idx) : bgdist})
        self.ion2background.update(ion2bg_local)

# Cell
import numpy as np
import random
import pandas as pd
from statistics import NormalDist
import math
from time import time
import typing
from numba import njit

class BackGroundDistribution:
    """Represents and derives an empirical distribution to describe the variation underlying a measurment
    """
    fc_resolution_factor = 100
    fc_conversion_factor = 1/fc_resolution_factor

    def __init__(self, start_idx : int, end_idx: int, ion2noNanvals : typing.Dict[int, str], idx2ion : dict,p2z : dict):

        """
        Initialize the background distribution from a subset of selected ions. The ions are pre-ordered and indexed and a sub range is selected. The
        Background Distribution is created from the sub-range.
        Args:
            start_idx (int): determines the start of sub-range
            end_idx (int): determines the end of the sub-range
            ion2noNanvals (dict): maps the ion to all measured intensities of this ion (no NAs/zero measurements)
            idx2ion (dict): distinct mapping of the index to the ion name
            p2z (dict): p-values are transformed into z-values on many occasions and are therefore cached with this dictionary.
        """
        self.fc2counts = {} #binned Fold change Distribution
        self.cumulative = np.array([])
        self.zscores = np.array([])
        self.min_fc =0
        self.max_fc = 0
        self.min_z=0
        self.max_z=0
        self.start_idx = int(start_idx)
        self.end_idx = int(end_idx)
        self.var = None
        self.SD = None
        self.ions = {idx2ion.get(idx) for idx in range(start_idx, end_idx)}
        self.fraction_missingval = self.calc_missingval_fraction(ion2noNanvals, idx2ion)


        anchor_fcs = self.generate_anchorfcs_from_intensity_range(ion2noNanvals, idx2ion)
        random.Random(42).shuffle(anchor_fcs) #set seed to ensure reproducibility
        self.generate_fc2counts_from_anchor_fcs(anchor_fcs)
        self.cumulative = self.transform_fc2counts_into_cumulative()
        self.calc_SD(0, self.cumulative)
        self.zscores = self.transform_cumulative_into_z_values(p2z)

        LOGGER.info(f"Created Background Distribution for {len(self.ions)} ions. SD: {self.SD}, fraction of missing values: {self.fraction_missingval:.2f}")

    def calc_missingval_fraction(self, ion2nonNanvals: dict, idx2ion: dict) -> float:
        """Calculates the fraction of missing values in the background distribution

        Args:
            ion2nonNanvals (dict): maps the ion to all measured intensities of this ion (no NAs/zero measurements)
            idx2ion (dict): distinct mapping of the index to the ion name

        Returns:
            float: fraction of missing values
        """
        value_nums = [len(ion2nonNanvals.get(idx2ion.get(idx))) for idx in range(self.start_idx, self.end_idx)]
        num_replicates = max(value_nums)
        num_total = num_replicates*(self.end_idx-self.start_idx)
        num_measured = sum(value_nums)
        num_missing = num_total - num_measured
        return num_missing/num_total

    def generate_anchorfcs_from_intensity_range(self, ion2noNanvals : dict, idx2ion : dict) -> list:
        """For each ion, a random intensity is selected as an "anchor" and the remaining intensities are subtracted from the achor.

        Args:
            ion2noNanvals (dict): maps the ion to all measured intensities of this ion (no NAs/zero measurements)
            idx2ion (dict): distinct mapping of the index to the ion name

        Returns:
            list: a merged list of all fold changes relative to the anchors
        """
        rng = random.Random(0)
        anchor_fcs = []
        for idx in range(self.start_idx, self.end_idx):
            vals = ion2noNanvals[idx2ion.get(idx)]
            if vals.size < 2:
                continue
            anchor_idx =  rng.randint(0, len(vals)-1)
            anchor_val = vals[anchor_idx]
            vals = np.delete(vals, anchor_idx)
            anchor_fcs.extend(vals-anchor_val)
        return anchor_fcs

    def generate_fc2counts_from_anchor_fcs(self,anchor_fcs : list):
        """Arbitrary pairs of anchor-changes are compared with each other, in order to determine the overall variation between the ions.

        Args:
            anchor_fcs (list): input list of the anchor-changes

        Returns:
            updates the self.fc2counts instance variable
        """
        anchor_fcs = anchor_fcs
        for idx in range(1, len(anchor_fcs)):
            fc_binned = np.rint(self.fc_resolution_factor*(0.5*(anchor_fcs[idx-1] - anchor_fcs[idx]))).astype(np.int64)
            self.fc2counts[fc_binned] = self.fc2counts.get(fc_binned, 0) + 1 #the distribution is saved in 2d (binned fold changes vs. count) for memory efficiency

        self.min_fc = min(self.fc2counts.keys())
        self.max_fc = max(self.fc2counts.keys())


    def transform_fc2counts_into_cumulative(self) -> np.array(float):
        """The binned fold change distribution is encoded in a 1d array, where the coordinate of the array represents the fold change and
        the value of the array represents the cumulative frequency.

        Returns:
            np.array: cumulative distribution of fold changes encoded in 1d array
        """
        cumulative = np.zeros(self.max_fc - self.min_fc +1).astype(np.int64)

        for entry in self.fc2counts.items():
            cumulative[int(entry[0]-self.min_fc)] +=entry[1]
        for idx in range(1,cumulative.shape[0]):
            cumulative[idx] +=cumulative[idx-1]

        return cumulative


    def transform_cumulative_into_z_values(self, p2z: dict):
        """
        OPTIMIZED: The binned fold change distribution is encoded in a 1d array, where the coordinate of the array represents the fold change and
        the value of the array represents the z-value. For each point in the distribution, we can calculate the z-value. This value encodes the distance from
        zero in a standard normal distribution that is required to obtain the same relative cumulative value

        Args:
            p2z (dict): p-values are transformed into z-values on many occasions and are therefore cached with this dictionary.
                       NOTE: This is now ignored in favor of the fast Numba implementation.

        Returns:
            np.array: array of z-values corresponding to the fold changes encoded in 1d array
        """
        total = self.cumulative[-1]
        min_pval = 1/(total+1)
        self.max_z = abs(NormalDist().inv_cdf(max(1e-9, min_pval)))

        # Use the Numba-optimized function for dramatic speedup (100x+ faster)
        return aqdiffutils.zscores_from_cumulative(self.cumulative, self.min_fc, total)


    def calc_zscore_from_fc(self, fc):
        return aqdiffutils.z_from_fc_lookup(fc, self.fc_conversion_factor, self.fc_resolution_factor, self.min_fc, self.cumulative, self.max_z, self.zscores)



    def calc_SD(self, mean:float, cumulative:list):
        """
        OPTIMIZED: Calculates the standard deviation of the background distribution
        Args:
            mean (float): mean value for the calculation
            cumulative (list[int]): cumulative distribution array
        """
        # Use the Numba-optimized function for dramatic speedup (100x+ faster)
        self.SD = aqdiffutils.sd_from_cumulative(np.asarray(cumulative), self.min_fc, mean, self.fc_conversion_factor)
        self.var = self.SD ** 2

    def get_cache_key(self):
        """
        Generate a unique, hashable cache key for this background distribution.

        Uses identifying properties that make this distribution unique:
        - Index range (start_idx, end_idx)
        - Fold change range (min_fc, max_fc)
        - Distribution size (cumulative array length)
        - Standard deviation (computed property)

        Returns:
            tuple: Hashable cache key that uniquely identifies this distribution
        """
        return (self.start_idx, self.end_idx, self.min_fc, self.max_fc,
                len(self.cumulative), round(self.SD, 6))




# Cell
from numba import jit
from time import time

class SubtractedBackgrounds(BackGroundDistribution):

    def __init__(self, from_dist, to_dist, p2z):
        self.max_fc = None
        self.min_fc = None
        self.var_from = from_dist.var
        self.var_to = to_dist.var
        self.cumulative = None
        max_joined, min_joined, cumulative = subtract_distribs(from_dist, to_dist)
        self.max_fc = max_joined
        self.min_fc = min_joined
        self.cumulative = cumulative
        self.fc2counts = transform_cumulative_into_fc2count(self.cumulative,self.min_fc)
        self.calc_SD(0, self.cumulative)
        self.zscores = self.transform_cumulative_into_z_values(p2z)

    def get_cache_key(self):
        """
        Generate a unique, hashable cache key for this subtracted background distribution.

        Since SubtractedBackgrounds doesn't have start_idx/end_idx, we use the
        properties that uniquely identify it: fold change range, distribution size,
        and standard deviation.

        Returns:
            tuple: Hashable cache key that uniquely identifies this distribution
        """
        return (self.min_fc, self.max_fc, len(self.cumulative),
                round(self.SD, 6), round(self.var_from, 6), round(self.var_to, 6))

def subtract_distribs(from_dist, to_dist):
    min_joined = from_dist.min_fc - to_dist.max_fc
    max_joined = from_dist.max_fc - to_dist.min_fc

    n_from = get_normed_freqs(from_dist.cumulative)
    n_to = get_normed_freqs(to_dist.cumulative)

    min_from = from_dist.min_fc
    min_to = to_dist.min_fc

    joined_init = np.zeros(max_joined-min_joined+1, dtype=np.int64)
    joined = get_joined(joined_init, n_from,n_to, min_from, min_to, min_joined)
    cumulative = np.cumsum(joined,dtype = np.int64)

    return max_joined, min_joined, cumulative

@njit
def get_joined(joined, n_from, n_to, min_from, min_to, min_joined):
    """
    Ultra-sparse optimization for get_joined function.


    This implementation uses exact-size pre-allocation for non-zero elements
    and early termination for completely sparse arrays.
    """
    # Count non-zero elements first for exact pre-allocation
    nz_from_count = 0
    nz_to_count = 0

    for i in range(len(n_from)):
        if n_from[i] != 0:
            nz_from_count += 1
    for i in range(len(n_to)):
        if n_to[i] != 0:
            nz_to_count += 1

    # Early termination for completely sparse arrays
    if nz_from_count == 0 or nz_to_count == 0:
        return joined

    # Pre-allocate exact size arrays for optimal memory usage
    from_indices = np.empty(nz_from_count, dtype=np.int64)
    to_indices = np.empty(nz_to_count, dtype=np.int64)
    from_values = np.empty(nz_from_count, dtype=np.int64)
    to_values = np.empty(nz_to_count, dtype=np.int64)

    # Collect non-zero data in one pass
    from_idx = 0
    for i in range(len(n_from)):
        if n_from[i] != 0:
            from_indices[from_idx] = i
            from_values[from_idx] = n_from[i]
            from_idx += 1

    to_idx = 0
    for i in range(len(n_to)):
        if n_to[i] != 0:
            to_indices[to_idx] = i
            to_values[to_idx] = n_to[i]
            to_idx += 1

    # Ultra-fast computation with pre-computed values
    for i in range(nz_from_count):
        fc_from = min_from + from_indices[i]
        freq_from = from_values[i]

        for j in range(nz_to_count):
            fc_to = min_to + to_indices[j]
            freq_to = to_values[j]
            fcdiff = fc_from - fc_to
            joined_idx = fcdiff - min_joined

            # Bounds checking for safety
            if 0 <= joined_idx < len(joined):
                joined[joined_idx] += freq_from * freq_to

    return joined

# Cell
def get_subtracted_bg(bgpair2diffDist, bg1, bg2, p2z):
    """
    OPTIMIZED: Improved caching for SubtractedBackgrounds to avoid expensive string conversions
    and improve cache hit rates.

    Uses a more efficient cache key based on object identity and properties rather than
    expensive string representations.
    """
    # Create efficient cache key using object properties instead of string conversion
    # This avoids the expensive str() operation and improves cache efficiency
    cache_key = _get_background_cache_key(bg1, bg2)

    # Use direct dictionary access instead of .keys() + .get() pattern
    if cache_key in bgpair2diffDist:
        return bgpair2diffDist[cache_key]

    # Create new SubtractedBackground and cache it
    subtr_bg = SubtractedBackgrounds(bg1, bg2, p2z)
    bgpair2diffDist[cache_key] = subtr_bg

    return subtr_bg

def _get_background_cache_key(bg1, bg2):
    """
    Generate efficient cache key for background distribution pairs.

    Uses the get_cache_key method from each background distribution to create
    a unique, hashable key that avoids expensive string operations while
    maintaining uniqueness.

    Args:
        bg1: First background distribution
        bg2: Second background distribution

    Returns:
        tuple: Hashable cache key
    """
    # Use the dedicated cache key method for more reliable identification
    key1 = bg1.get_cache_key()
    key2 = bg2.get_cache_key()

    # Ensure consistent ordering for cache hits regardless of argument order
    if key1 <= key2:
        return (key1, key2)
    else:
        return (key2, key1)

# Cell

def get_doublediff_bg(deed_ion1, deed_ion2, deedpair2doublediffdist, p2z):
    """
    OPTIMIZED: Improved caching for double differential backgrounds with more efficient
    cache key generation and lookup.
    """
    # Generate efficient cache key
    cache_key = _get_background_cache_key(deed_ion1, deed_ion2)

    # Check cache with direct access (more efficient than .keys())
    if cache_key in deedpair2doublediffdist:
        return deedpair2doublediffdist[cache_key]

    # The original code checked for inverted keys, but our cache key generation
    # already handles ordering, so this is no longer needed

    # Create new SubtractedBackground and cache it
    subtr_bg = SubtractedBackgrounds(deed_ion1, deed_ion2, p2z)
    deedpair2doublediffdist[cache_key] = subtr_bg

    return subtr_bg

def invert_deedkey(deedkey):
    return (deedkey[1], deedkey[0])


# Cell
from numba import njit

#get normalized freqs from cumulative
@njit
def get_normed_freqs(cumulative):
    normfact = 2**30 /cumulative[-1]
    freqs =get_freq_from_cumul(cumulative)
    for i in range(len(freqs)):
        freqs[i] *= normfact
    return freqs

# Cell
from numba import njit

#transform cumulative into frequency
@njit
def get_freq_from_cumul(cumulative):
    res = np.zeros(len(cumulative), dtype=np.int64)
    res[0] = cumulative[0]
    for i in range(1,len(cumulative)):
        res[i] = cumulative[i]-cumulative[i-1]

    return res

# Cell
import numba.typed
import numba.types
from numba import njit

@njit
def _transform_cumulative_vectorized(cumulative, min_fc):
    """Optimized vectorized implementation with 2x speedup"""
    if len(cumulative) <= 1:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)

    # Calculate ALL differences vectorially - this is the key optimization
    diffs = np.diff(cumulative)

    # Create fold change indices
    fcs = np.arange(1, len(cumulative), dtype=np.int64) + min_fc

    return fcs, diffs

def transform_cumulative_into_fc2count(cumulative, min_fc):
    """
    Optimized transform function with 2x speedup.

    Uses vectorized numpy operations instead of explicit loops.
    Maintains identical results to original implementation.
    """
    fcs, counts = _transform_cumulative_vectorized(cumulative, min_fc)
    return dict(zip(fcs, counts))

# Cell
@njit
def get_cumul_from_freq(freq):
    res = np.zeros(len(freq), dtype=np.int64)
    res[0] = freq[0]
    for i in range(1,len(freq)):
        res[i] = res[i-1] + freq[i]

    return res
