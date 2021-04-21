import alphaquant.diffquant_utils as aqutils

tabledir = "../../test_data/input_table_formats/"
output_folder = "./results"
import pkg_resources

print(pkg_resources.resource_filename("alphaquant", "core.py"))
import pathlib
print(pathlib.Path().absolute())

aqutils.import_data(f"{tabledir}/diann_precursor.tsv", samplemap_file=f"{tabledir}/samplemap.diann.tsv",output_folder=output_folder)