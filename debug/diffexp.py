# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
# default_exp background_distributions


# %%
class ConditionBackgrounds():

    def __init__(self, normed_condition_df):
        self.backgrounds = []
        self.normed_condition_df = normed_condition_df
        self.ion2background = {}
        self.ion2nonNanvals = {}
        self.init_ion2nonNanvals(self.normed_condition_df)
        self.select_intensity_ranges()

    def init_ion2nonNanvals(self, normed_condition_df):
        for peptide, vals in normed_condition_df.iterrows():
            self.ion2nonNanvals[peptide] = vals[vals.notna()].values


    def select_intensity_ranges(self):
        total_available_comparisons =0
        num_contexts = 100
        cumulative_counts = np.zeros(self.normed_condition_df.shape[0])

        for idx ,count in enumerate(self.normed_condition_df.count(axis=1)):
            total_available_comparisons+=count-1
            cumulative_counts[idx] = int(total_available_comparisons/2)
        
        
        #assign the context sizes
        context_size = np.max([1000, int(total_available_comparisons/(1+num_contexts/2))])
        halfcontext_size = int(context_size/2)
        context_boundaries = np.zeros(3).astype(int)

        middle_idx = int(np.searchsorted(cumulative_counts, halfcontext_size))
        end_idx = int(np.searchsorted(cumulative_counts, context_size))


        context_boundaries[0] = 0
        context_boundaries[1] = middle_idx
        context_boundaries[2] = end_idx
        while context_boundaries[1] < len(cumulative_counts)-1:
            self.backgrounds.append(BackGroundDistribution(context_boundaries[0], context_boundaries[2], self.ion2nonNanvals))
            context_boundaries[0] = context_boundaries[1]
            context_boundaries[1] = context_boundaries[2]
            end_idx = np.searchsorted(cumulative_counts, context_size + cumulative_counts[context_boundaries[0]])
            if end_idx > len(cumulative_counts)-(context_boundaries[1]-context_boundaries[0])/1.5:
                end_idx = len(cumulative_counts)-1
            context_boundaries[2] = end_idx




# %%
import numpy as np
from random import shuffle
import pandas as pd
from scipy.stats import norm
import math

class BackGroundDistribution:

    fc_resolution_factor = 100
    fc_conversion_factor = 1/fc_resolution_factor

    def __init__(self, start_idx, end_idx, ion2noNanvals):
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

        anchor_fcs = self.generate_anchorfcs_from_intensity_range(ion2noNanvals)
        shuffle(anchor_fcs)
        self.generate_fc2counts_from_anchor_fcs(anchor_fcs)
        self.cumulative = self.transform_fc2counts_into_cumulative()
        self.calc_SD(0, self.cumulative)
        self.zscores = self.transform_cumulative_into_z_values()

    def generate_anchorfcs_from_intensity_range(self,ion2noNanvals):

        anchor_fcs = []
        for idx in range(self.start_idx, self.end_idx):
            vals = ion2noNanvals[idx]
            if vals.size < 2:
                continue
            anchor_idx =  np.random.randint(0, len(vals))
            anchor_val = vals[anchor_idx]
            vals = np.delete(vals, anchor_idx)
            anchor_fcs.extend(vals-anchor_val)
        return anchor_fcs


    def generate_fc2counts_from_anchor_fcs(self,anchor_fcs):
        
        anchor_fcs = np.array(anchor_fcs)
        for idx in range(1, anchor_fcs.shape[0]):
            fc_binned = np.rint(self.fc_resolution_factor*(0.5*(anchor_fcs[idx-1] - anchor_fcs[idx]))).astype(np.long)
            self.fc2counts[fc_binned] = self.fc2counts.setdefault(fc_binned, 0) + 1

        self.min_fc = min(self.fc2counts.keys())
        self.max_fc = max(self.fc2counts.keys())

    
    def transform_fc2counts_into_cumulative(self):
        
        cumulative = np.zeros(self.max_fc - self.min_fc +1).astype(np.long)

        for entry in self.fc2counts.items():
            cumulative[int(entry[0]-self.min_fc)] +=entry[1]
        for idx in range(1,cumulative.shape[0]):
            cumulative[idx] +=cumulative[idx-1]
        
        return cumulative

    
    def transform_cumulative_into_z_values(self):
        total = self.cumulative[-1]
        min_pval = 1/(total+1)
        self.max_z = abs(norm.ppf(max(1e-9, min_pval)))
        zscores = np.empty(len(self.cumulative))
        zero_pos = -self.min_fc

        normfact_posvals = 1/(total-self.cumulative[zero_pos]+1)
        normfact_negvals = 1/(self.cumulative[zero_pos-1]+1)
        for i in range(len(self.cumulative)):
            num_more_extreme = 0
            if i == zero_pos:
                zscores[i] = 0
                continue
            if i!=zero_pos and i<len(self.cumulative)-1:
                num_more_extreme = self.cumulative[i] if i<zero_pos else  self.cumulative[-1] - self.cumulative[i+1]
            
            normfact = normfact_negvals if i<zero_pos else normfact_posvals
            p_val = 0.5*max(1e-9, (num_more_extreme+1)*normfact)
            sign = -1 if i<zero_pos else 1
            zscores[i] = sign*norm.ppf(p_val) ##ppf is the inverese cumulative distribution function

        return zscores


    def calc_zscore_from_fc(self, fc):
        if abs(fc)<1e-9:
            return 0
        k = int(fc * self.fc_resolution_factor)
        print(f"type k is {type(k)}")
        rank = k-self.min_fc
        print(f"minfc type is {type(self.min_fc)}")
        if rank <0:
            return -self.max_z
        if rank >=len(self.cumulative):
            return self.max_z
        print(rank)
        return self.zscores[rank]


    def calc_SD(self, mean, cumulative):
        sq_err = 0.0
        previous =0
        for i in range(len(cumulative)):
            fc = (i+self.min_fc)*self.fc_conversion_factor
            sq_err += (cumulative[i] - previous)*(fc-mean)**2
            previous = cumulative[i]
        total = cumulative[-1]
        var = sq_err/total
        self.var = var
        print(f"var is {var}")
        self.SD = math.sqrt(var)


# %%
#hide
#test background distribution
from scipy.stats import norm
import matplotlib.pyplot as plt

idx2nonnanvals = {}

for idx in range(100000):
    nonnanvals =  np.random.normal(loc=0, size=3)
    idx2nonnanvals[idx] = nonnanvals
    
bgdist = BackGroundDistribution(0, 99999, idx2nonnanvals)

def tranform_fc2count_to_fc_space(fc2counts, num_fcs, rescale_factor):
    fc2counts_fcscales = {}
    for fc, count in fc2counts.items():
        fc2counts_fcscales[fc*rescale_factor] = count/num_fcs

    return fc2counts_fcscales

fc2counts_rescaled = tranform_fc2count_to_fc_space(bgdist.fc2counts, bgdist.cumulative[-1],1/100.0)

plt.bar(list(fc2counts_rescaled.keys()), fc2counts_rescaled.values(),width=0.01,color='g',)
axes2 = plt.twinx()
x = np.linspace(-4, 4, 1000)
axes2.plot(x, norm.pdf(x, 0, 1)/1.15)
axes2.set_ylim(0.0, 0.4)
plt.show()



# %%
#subtract two Empirical Backgrounds
from scipy.stats import norm
class SubtractedBackgrounds(BackGroundDistribution):

    def __init__(self, from_dist, to_dist):
        self.max_fc = None
        self.min_fc = None
        self.cumulative = None
        self.subtract_distribs(from_dist, to_dist)
        self.calc_SD(0, self.cumulative)
        self.zscores = self.transform_cumulative_into_z_values()
        
    def subtract_distribs(self,from_dist, to_dist):
        min_joined = from_dist.min_fc - to_dist.max_fc
        max_joined = from_dist.max_fc - to_dist.min_fc

        n_from = get_normed_freqs(from_dist.cumulative)
        n_to = get_normed_freqs(to_dist.cumulative)

        min_from = from_dist.min_fc
        min_to = to_dist.min_fc

        joined = np.empty(max_joined-min_joined+1, dtype="long")
        
        for from_idx in range(len(n_from)):
            fc_from = min_from + from_idx
            freq_from = n_from[from_idx]
            for to_idx in range(len(n_to)):
                fc_to = min_to + to_idx
                freq_to = n_to[to_idx]
                fcdiff = fc_from - fc_to
                joined_idx = fcdiff - min_joined
                if freq_from*freq_to < 0:
                    print("smaller zero")
                joined[joined_idx] += (freq_from*freq_to)
        self.max_fc = max_joined
        self.min_fc = min_joined
        self.cumulative = joined
 


# %%
#evaluate differential regulation peptide
from scipy.stats import norm
import numpy as np
import math

class DifferentialIon():
    
    def __init__(self,noNanvals_from, mean_from, bgdist_from, noNanvals_to, mean_to,  bgdist_to):
        self.usable = False
        p_val, fc, z_val = calc_diffreg_peptide(noNanvals_from, mean_from,bgdist_from, noNanvals_to, mean_to, bgdist_to)
        if (p_val!=None):
            self.p_val=p_val
            self.fc=fc
            self.z_val = z_val
            self.usable = True

def calc_diffreg_peptide(noNanvals_from, mean_from,bgdist_from, noNanvals_to, mean_to,  bgdist_to): #TO Do normalize the input vectors between conds
    nrep_from = len(noNanvals_from)
    nrep_to = len(noNanvals_to)

    if ((nrep_from==0) or (nrep_to ==0)):
        return None, None, None

    diffDist = SubtractedBackgrounds(bgdist_from, bgdist_to)

    perEvidenceVariance = diffDist.var + (nrep_to-1) * bgdist_from.var + (nrep_from-1 ) * bgdist_to.var
    totalVariance = perEvidenceVariance*nrep_to * nrep_from
    fc_sum =0
    z_sum=0

    for from_intens in noNanvals_from:
        for to_intens in noNanvals_to:
            fc = from_intens - to_intens
            fc_sum+=fc
            z_sum += diffDist.calc_zscore_from_fc(fc)

    fc = fc_sum/(nrep_from * nrep_to)
    outlier_scaling_factor = calc_outlier_scaling_factor(noNanvals_from, mean_from, noNanvals_to, mean_to, bgdist_from, bgdist_to, diffDist)
    scaled_SD = outlier_scaling_factor * math.sqrt(totalVariance/diffDist.var)
    p_val = 2.0 * (1.0 -  norm(loc=0, scale= scaled_SD).cdf(abs(z_sum)))
    z_val = z_sum/scaled_SD
    return p_val, fc, z_val

def calc_outlier_scaling_factor(noNanvals_from, mean_from, noNanvals_to, mean_to, bgdist_from, bgdist_to, diffDist):
    between_rep_SD_from = math.sqrt(sum(np.square(noNanvals_from-mean_from))) if len(noNanvals_from)>1 else bgdist_from.SD
    between_rep_SD_to = math.sqrt(sum(np.square(noNanvals_to-mean_to))) if len(noNanvals_to)>1 else bgdist_to.SD

    highest_SD_from = max(between_rep_SD_from, bgdist_from.SD)
    highest_SD_to = max(between_rep_SD_to, bgdist_to.SD)
    highest_SD_combined = math.sqrt(highest_SD_from**2 + highest_SD_to**2)

    scaling_factor = max(1.0, highest_SD_combined/diffDist.SD)
    return scaling_factor



# %%
num_vals = 10000
rand1 = np.random.normal(loc=0, size=num_vals)
rand2 = np.random.normal(loc=0, size=num_vals)
rand3 = np.random.normal(loc=0, size=num_vals)
rand4 = np.random.normal(loc=0, size=num_vals)
rand5 = np.random.normal(loc=0, size=num_vals)

randarray = pd.DataFrame({1:rand1, 2:rand2, 3:rand3, 4:rand4, 5:rand5})
#display(randarray)
condbg = ConditionBackgrounds(randarray)


# %%
#test

def test_subtract_distribs():
    from_dist = [1,1,2,1,1]
    to_dist = [1,1,2,1,1]


# %%
#hide
#test subtract background distribution
from scipy.stats import norm
import matplotlib.pyplot as plt

idx2nonnanvals = {}

for idx in range(20000):
    nonnanvals =  np.random.normal(loc=0, size=3)
    idx2nonnanvals[idx] = nonnanvals
    
bgdist1 = BackGroundDistribution(0, 9999, idx2nonnanvals)
bgdist2 = BackGroundDistribution(10000, 19999, idx2nonnanvals)

subtracted_bgs = SubtractedBackgrounds(bgdist1, bgdist2)

def tranform_fc2count_to_fc_space(fc2counts, num_fcs, rescale_factor):
    fc2counts_fcscales = {}
    for fc, count in fc2counts.items():
        fc2counts_fcscales[fc*rescale_factor] = count/num_fcs

    return fc2counts_fcscales

fc2counts_rescaled = tranform_fc2count_to_fc_space(subtracted_bgs.fc2counts, bgdist.cumulative[-1],1/100.0)

plt.bar(list(fc2counts_rescaled.keys()), fc2counts_rescaled.values(),width=0.01,color='g',)
axes2 = plt.twinx()
x = np.linspace(-4, 4, 1000)
axes2.plot(x, norm.pdf(x, 0, 2)/1.15)
axes2.set_ylim(0.0, 0.4)
plt.show()


# %%
import math
import statistics
from scipy.stats import norm

class DifferentialProtein():

    def __init__(self, name):
        self.pval=None
        self.fc=None
        self.name = name

    def evaluate_protein_expression(self, ion_diffresults):
        ion_diffresults = list(filter(lambda _f : _f.usable, ion_diffresults))
        if len(ion_diffresults) ==0:
            return
        fcs = list(map(lambda _dr : _dr.fc,ion_diffresults))
        median_fc = statistics.median(fcs)
        ion_diffresults = select_robust_if_many_ions(fcs, median_fc,ion_diffresults)
        z_sum = sum(map(lambda _dr: _dr.zval, table))
        p_val = 2.0 * (1.0 - norm(0, math.sqrt(len(ion_diffresults))).cdf(abs(z_sum)))

        return median_fc, p_val


    def select_robust_if_many_ions(self, fcs, median_fc,ion_diffresults):
        ion_diffresults = sort(lambda _dr : abs(_dr.fc - median_fc),ion_diffresults)
        ninety_perc_cutoff = math.ceil(0.9*len(ion_diffresults)) #the ceil function ensures that ions are only excluded if there are more than 10 available
        if ninety_perc_cutoff >0:
            ion_diffresults = ion_diffresults[:ninety_perc_cutoff]
        return ion_diffresults


# %%
#transform cumulative into frequency

def get_freq_from_cumul(cumulative):
    res = np.empty(len(cumulative), dtype="long")
    res[0] = cumulative[0]
    for i in range(1,len(cumulative)):
        res[i] = cumulative[i]-cumulative[i-1]

    return res


# %%
#test
def test_get_freq_from_cumul():
            arr = [1,2,3,4,5,6]
            freqs = get_freq_from_cumul(arr)
            assert (freqs == [1,1,1,1,1,1]).all()
test_get_freq_from_cumul()


# %%
#get normalized freqs from cumulative

def get_normed_freqs(cumulative):
    normfact = 2**30 /cumulative[len(cumulative)-1]
    freqs =get_freq_from_cumul(cumulative)
    for i in range(len(freqs)):
        freqs[i] *= normfact
    return freqs


# %%
#test
def test_diffreg_pep(condbg):
    noNanvals_from = np.array([1,2,1,2])
    noNanvals_to = np.array([2,1,2,1])
    mean_from =1.5
    mean_to = 1.5
    bgdist_from = condbg.backgrounds[0]
    bgdist_to = condbg.backgrounds[1]
    diffion = DifferentialIon(noNanvals_from,mean_from, bgdist_from, noNanvals_to,mean_to, bgdist_to)
    print(f"fc {diffion.fc}, pval {diffion.p_val}, zval {diffion.z_val} ")
test_diffreg_pep(condbg)


# %%



