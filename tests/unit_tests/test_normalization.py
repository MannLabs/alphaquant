import alphaquant.norm.normalization as aq_norm
import numpy as np

def test_merged_distribs():
    anchor_distrib = np.array([1, 1, 1, 1, 1])
    shift_distrib = np.array([2, 2, 2, 2, 2])
    counts_anchor_distrib = 4
    counts_shifted_distib = 1
    assert (aq_norm.merge_distribs(anchor_distrib, shift_distrib, counts_anchor_distrib, counts_shifted_distib)== np.array([1.2, 1.2, 1.2, 1.2, 1.2])).any()
    print("test_merged_distribs passed")

test_merged_distribs() 


def generate_randarrays(number_arrays,size_of_array):
    randarray = []
    for i in range(number_arrays):
        shift = np.random.uniform(low=-10, high=+10)
        randarray.append(np.random.normal(loc=shift, size=size_of_array))
    return np.array(randarray)

import numpy as np

def test_sampleshift():
    randarray = generate_randarrays(5, 1000)
    sample2shift = aq_norm.get_normfacts_withincond(randarray)
    normalized_randarray = aq_norm.apply_sampleshifts(randarray, sample2shift)
    num_samples = normalized_randarray.shape[0]
    merged_sample = []
    for i in range(num_samples):
        merged_sample.extend(normalized_randarray[i])
    stdev = np.std(merged_sample)
    print(f"STDev {stdev}")
    assert (stdev <=1.2) 
    


test_sampleshift()