import sys
sys.path.append('/Users/constantin/workspace/MS-EmpiRe_Python/')
from ms_empire.background_distributions import *
from ms_empire.normalization import *
from ms_empire.diff_analysis import *
from ms_empire.visualizations import *
from ms_empire.benchmarking import *
from ms_empire.diffquant_utils import *
from ms_empire.diff_analysis_manager import *

import os
os.chdir("/Users/constantin/workspace/MS-EmpiRe_Python")

import os 
import matplotlib.pyplot as plt

peptides_tsv = "./test_data/ap_quant_peptides.tsv"
samplemap_tsv = "./test_data/ap.samples.map"
unnormed_df, labelmap_df = read_tables(peptides_tsv, samplemap_tsv, "precursor")
unnormed_df = unnormed_df.set_index("protein")
protvals = unnormed_df.loc["sp|Q9Y2R4|DDX52_HUMAN"].dropna(thresh=1, axis=0)
protvals = protvals.to_numpy()
print(protvals)
normed_protvals = normalize_withincond(protvals)
plt.plot(protvals,label = 'protvals')
plt.show()
plt.plot(normed_protvals,label = 'protvals')
display(normed_df)