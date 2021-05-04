import alphaquant.diffquant_utils as aqutils
import alphaquant.diff_analysis_manager as diffmgr
import os
os.chdir("/Users/constantin/workspace/EmpiRe/alphaquant/tests/input_parsing")
tabledir = os.path.join("..", "..", "test_data", "input_table_formats")
results_folder = os.path.join(".", "results")
print(os.path.abspath(results_folder))


#input_file = os.path.join(tabledir,"diann.tsv" )
input_file = os.path.join(tabledir, "spectronaut.tsv")
#input_file = os.path.join(tabledir, "mq_peptides.txt")

#samplemap_file = os.path.join(tabledir, "samplemap.diann.tsv")
samplemap_file = os.path.join(tabledir, "samplemap.spectronaut.tsv")
#samplemap_file = os.path.join(tabledir, "samplemap.mq.tsv")

#import the input table once the input and the results folder are specified. 
# The function automatically recognizes the format (Currently MQ, Spectronaut, DIA-NN configured)
input_data = aqutils.import_data(input_file,results_folder=results_folder)

#get sample names from the imported table
samplenames = aqutils.get_samplenames(input_data)

#load the samplemap dataframe (in case the user uploads a file. Basically a pandas import + separator check)
samplemap_df = aqutils.load_samplemap(samplemap_file)

#compare samplemap and actual table, merge & logtransform intensities
input_processed, samplemap_df_processed = aqutils.prepare_loaded_tables(input_data, samplemap_df)



#run the pipeline, plots are deactivated (runtime_plots actually set False as default also)
diffmgr.run_pipeline(input_processed, samplemap_df_processed, results_folder,runtime_plots=False)