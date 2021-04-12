import sys
sys.path.append('/Users/constantin/workspace/MS-EmpiRe_Python/')
from alphaquant.background_distributions import *
from alphaquant.normalization import *
from alphaquant.diff_analysis import *
from alphaquant.visualizations import *
from alphaquant.benchmarking import *
from alphaquant.diffquant_utils import *
from alphaquant.diff_analysis_manager import *
from alphaquant.protein_intensities import *

import os
os.chdir("/Users/constantin/workspace/MS-EmpiRe_Python")



import pandas as pd

mq_lfq = pd.read_csv("./test_data/proteinGroups.txt", sep = "\t", usecols = ["LFQ intensity Shotgun_12-01_1", "LFQ intensity Shotgun_12-01_2", "LFQ intensity Shotgun_12-01_3",
"LFQ intensity Shotgun_02-01_1", "LFQ intensity Shotgun_02-01_2", "LFQ intensity Shotgun_02-01_3", "Species", "id"]).rename(columns = {"Species":"species",
"id" : "protein"})
mq_lfq["species"] = mq_lfq["species"].replace(to_replace={"Homo sapiens" : "Human", "Escherichia coli" : "Ecoli"})
mq_lfq = mq_lfq.astype({"protein" :"str"})
minrep = 1

ap_c1_filt_idx = mq_lfq[["LFQ intensity Shotgun_12-01_1", "LFQ intensity Shotgun_12-01_2", "LFQ intensity Shotgun_12-01_3"]].dropna(thresh=minrep, axis=0).index
ap_c2_filt_idx = mq_lfq[["LFQ intensity Shotgun_02-01_1", "LFQ intensity Shotgun_02-01_2", "LFQ intensity Shotgun_02-01_3"]].dropna(thresh=minrep, axis=0).index

mq_lfq["HeLa12_median_ref"] = mq_lfq[["LFQ intensity Shotgun_12-01_1", "LFQ intensity Shotgun_12-01_2", "LFQ intensity Shotgun_12-01_3"]].median(axis = 1)
mq_lfq["HeLa2_median_ref"] = mq_lfq[["LFQ intensity Shotgun_02-01_1", "LFQ intensity Shotgun_02-01_2", "LFQ intensity Shotgun_02-01_3"]].median(axis = 1)
idx_intersect = ap_c1_filt_idx.intersection(ap_c2_filt_idx)
mq_lfq = mq_lfq.loc[idx_intersect]
plot_fold_change(mq_lfq, "LFQ intensity Shotgun_02-01_1", "LFQ intensity Shotgun_12-01_2")
