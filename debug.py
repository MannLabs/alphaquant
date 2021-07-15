import alphaquant.diff_analysis_manager as aqmgr
import alphaquant.diffquant_utils as aqutils

# fragion_df = aqutils.import_data("/Users/constantin/workspace/Quantification/Benchmarks/DIA/10_1_and_100_1/small_subset/benchmark_subset.tsv.aq_reformat.tsv")
# samplemap = aqutils.load_samplemap("/Users/constantin/workspace/Quantification/Benchmarks/DIA/10_1_and_100_1/small_subset/samples.map.tsv")
# condpair_combinations = [('Y1', 'Y10'), ('Y1', 'Y100')]
# fragion_df, samplemap = aqutils.prepare_loaded_tables(fragion_df, samplemap)
# aqmgr.run_pipeline(fragion_df, samplemap, condpair_combinations=condpair_combinations, runtime_plots=True)

# fragion_df = aqutils.import_data("/Users/constantin/workspace/Quantification/Benchmarks/MaxQuant/peptides.txt.aq_reformat.tsv")
# samplemap = aqutils.load_samplemap("/Users/constantin/workspace/Quantification/Benchmarks/MaxQuant/samples.map.tsv")
# fragion_df, samplemap = aqutils.prepare_loaded_tables(fragion_df, samplemap)

# aqmgr.run_pipeline(fragion_df, samplemap, runtime_plots=True)

#whole dataset
import alphaquant.diff_analysis_manager as aqmgr
import alphaquant.diffquant_utils as aqutils
 
fragion_df = aqutils.import_data("/Users/constantin/workspace/Quantification/Benchmarks/DIA/10_1_and_100_1/eval_noise/yeast_report.tsv.aq_reformat.tsv")
samplemap = aqutils.load_samplemap("/Users/constantin/workspace/Quantification/Benchmarks/DIA/10_1_and_100_1/eval_noise/samples.map.swapped.tsv")
condpair_combinations = [('Y1', 'Y10'), ('Y1', 'Y100')]
fragion_df, samplemap = aqutils.prepare_loaded_tables(fragion_df, samplemap)
aqmgr.run_pipeline(fragion_df, samplemap, condpair_combinations=condpair_combinations, minrep = 5, runtime_plots=True)