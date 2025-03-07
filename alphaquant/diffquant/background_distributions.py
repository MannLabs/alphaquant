
from time import time
import numpy as np
import alphaquant.diffquant.diffutils as aqutils

import alphaquant.config.config as aqconfig
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)


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
        t_start = time()
        self.ion2nonNanvals = aqutils.get_non_nas_from_pd_df(normed_condition_df)
        self.ion2allvals = aqutils.get_ionints_from_pd_df(normed_condition_df)
        t_end = time()
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


    def transform_cumulative_into_z_values(self:int, p2z: dict):
        """
        The binned fold change distribution is encoded in a 1d array, where the coordinate of the array represents the fold change and
        the value of the array represents the z-value. For each point in the distribution, we can calculate the z-value. This value encodes the distance from
        zero in a standard normal distribution that is required to obtain the same relative cumulative value

        Args:
            p2z (dict): p-values are transformed into z-values on many occasions and are therefore cached with this dictionary.

        Returns:
            np.array: array of z-values corresponding to the fold changes encoded in 1d array
        """
        total = self.cumulative[-1]
        min_pval = 1/(total+1)
        self.max_z = abs(NormalDist().inv_cdf(max(1e-9, min_pval)))
        zscores = np.zeros(len(self.cumulative))
        zero_pos = -self.min_fc

        normfact_posvals = 1/(total-self.cumulative[zero_pos]+1)
        normfact_negvals = 1/(self.cumulative[zero_pos-1]+1)
        for i in range(len(self.cumulative)):
            t_start = time()
            num_more_extreme = 0
            normfact = 0
            if i == zero_pos or i==len(self.cumulative)-1:
                zscores[i] = 0
                continue

            if i < zero_pos:
                num_more_extreme = self.cumulative[i]
                normfact = normfact_negvals
            else:
                num_more_extreme = self.cumulative[-1] - self.cumulative[i+1]
                normfact = normfact_posvals

            p_val = 0.5*max(1e-9, (num_more_extreme+1)*normfact)
            sign = -1 if i<zero_pos else 1
            t_empirical = time()
            zscore = sign*abs(get_z_from_p_empirical(p_val, p2z))
            zscores[i] =  zscore
            t_nd_lookup = time()
        return zscores


    def calc_zscore_from_fc(self, fc):
        return _calc_zscore_from_fc(fc, self.fc_conversion_factor, self.fc_resolution_factor, self.min_fc, self.cumulative, self.max_z, self.zscores)



    def calc_SD(self, mean:float, cumulative:list):
        """
        Calculates the standard deviation of the background distribution
        Args:
            mean (float): [description]
            cumulative (list[int]): [description]
        """
        sq_err = 0.0
        previous =0
        for i in range(len(cumulative)):
            fc = (i+self.min_fc)*self.fc_conversion_factor
            sq_err += (cumulative[i] - previous)*(fc-mean)**2
            previous = cumulative[i]
        total = cumulative[-1]
        var = sq_err/total
        self.var = var
        self.SD = math.sqrt(var)

@njit
def _calc_zscore_from_fc(fc, fc_conversion_factor, fc_resolution_factor, min_fc, cumulative, max_z, zscores):
    """
    Quick conversion function that looks up the z-value corresponding to an observed new fold change.
    The fold change is mapped to its fc-bin in the binned fold change distribution and then the z-value of the bin is looked up

    Args:
        fc (float): [description]

    Returns:
        float: z-value of the observed fold change, based on the background distribution
    """
    if abs(fc)<fc_conversion_factor:
        return 0
    k = int(fc * fc_resolution_factor)
    rank = k-min_fc
    if rank <0:
        return -max_z
    if rank >=len(cumulative):
        return max_z
    return zscores[rank]


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
        t_start = time()
        self.fc2counts = transform_cumulative_into_fc2count(self.cumulative,self.min_fc)
        t_cumul_transf = time()
        self.calc_SD(0, self.cumulative)
        t_calc_SD = time()
        self.zscores = self.transform_cumulative_into_z_values(p2z)
        t_calc_zvals = time()

def subtract_distribs(from_dist, to_dist):
    min_joined = from_dist.min_fc - to_dist.max_fc
    max_joined = from_dist.max_fc - to_dist.min_fc

    n_from = get_normed_freqs(from_dist.cumulative)
    n_to = get_normed_freqs(to_dist.cumulative)

    min_from = from_dist.min_fc
    min_to = to_dist.min_fc

    joined_init = np.zeros(max_joined-min_joined+1, dtype=np.int64)
    t_start = time()
    joined = get_joined(joined_init, n_from,n_to, min_from, min_to, min_joined)
    t_join = time()
    cumulative = np.cumsum(joined,dtype = np.int64)
    t_cumul = time()

    return max_joined, min_joined, cumulative

@jit(nopython=True)
def get_joined(joined,n_from, n_to, min_from, min_to, min_joined):
    count_comparisons =0
    for from_idx in range(len(n_from)):
        fc_from = min_from + from_idx
        freq_from = n_from[from_idx]
        for to_idx in range(len(n_to)):
            fc_to = min_to + to_idx
            freq_to = n_to[to_idx]
            fcdiff = fc_from - fc_to
            joined_idx = fcdiff - min_joined
            freq_multiplied = freq_from*freq_to
            joined[joined_idx] += (freq_multiplied)
            count_comparisons+=1
    return joined

# Cell
def get_subtracted_bg(bgpair2diffDist, bg1, bg2, p2z):

    bgpair = (str(bg1), str(bg2))
    if bgpair in bgpair2diffDist.keys():
        return bgpair2diffDist.get(bgpair)


    subtr_bg = SubtractedBackgrounds(bg1, bg2, p2z)
    bgpair2diffDist[bgpair] = subtr_bg

    return subtr_bg

# Cell

def get_doublediff_bg(deed_ion1, deed_ion2, deedpair2doublediffdist, p2z):

    deedkey = (str(deed_ion1), str(deed_ion2))
    inverted_deedkey = invert_deedkey(deedkey)

    if deedkey in deedpair2doublediffdist.keys():
        return deedpair2doublediffdist.get(deedkey)

    if inverted_deedkey in deedpair2doublediffdist.keys():
        return deedpair2doublediffdist.get(inverted_deedkey)

    subtr_bg = SubtractedBackgrounds(deed_ion1, deed_ion2, p2z)
    deedpair2doublediffdist[deedkey] = subtr_bg

    return subtr_bg

def invert_deedkey(deedkey):
    return (deedkey[1], deedkey[0])

# Cell
from statistics import NormalDist

def get_z_from_p_empirical(p_emp,p2z):
    p_rounded = np.format_float_scientific(p_emp, 1)
    if p_rounded in p2z:
        return p2z.get(p_rounded)
    z = NormalDist().inv_cdf(float(p_rounded))
    p2z[p_rounded] = z
    return z

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
#@njit

def transform_cumulative_into_fc2count(cumulative, min_fc):
#     res_dict = numba.typed.Dict.empty(
#     key_type=numba.types.int64,
#     value_type=numba.types.int64,
# )
    res_dict = {}
    for idx in range(1, len(cumulative)):
        fc = idx + min_fc
        res_dict[fc] = cumulative[idx] - cumulative[idx-1]
    return res_dict

# Cell
@njit
def get_cumul_from_freq(freq):
    res = np.zeros(len(freq), dtype=np.int64)
    res[0] = freq[0]
    for i in range(1,len(freq)):
        res[i] = res[i-1] + freq[i]

    return res
