import alphaquant.diffquant_utils as aqutils
import os
tabledir = "/Users/constantin/workspace/EmpiRe/alphaquant/test_data/input_table_formats"
results_folder = "/Users/constantin/workspace/EmpiRe/alphaquant/tests/input_parsing"
import pkg_resources

print(pkg_resources.resource_filename("alphaquant", "core.py"))
import pathlib
print(pathlib.Path().absolute())

#aqutils.import_data(f"{tabledir}/diann.tsv", samplemap_file=f"{tabledir}/samplemap.diann.tsv",results_folder=results_folder)

#aqutils.import_data(f"{tabledir}/spectronaut.tsv", samplemap_file=f"{tabledir}/samplemap.spectronaut.tsv",results_folder=results_folder)

aqutils.import_data(f"{tabledir}/mq_peptides.txt", samplemap_file=f"{tabledir}/samplemap.spectronaut.tsv",results_folder=results_folder)