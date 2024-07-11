import sys
import alphaquant.benchm.testfile_handling
test_folder = "../test_data/"
links_yaml_all_testfiles = "../alphaquant/config/download_links_for_testfiles_all.yaml"
links_yaml_quicktest_files = "../alphaquant/config/download_links_for_testfiles_quicktest.yaml"



if __name__ == '__main__':
    command_line_arguments = sys.argv
    type_of_test = command_line_arguments[1]

    if type_of_test == 'quicktest':
        links_yaml = links_yaml_quicktest_files
    elif type_of_test == 'all_tests':
        links_yaml = links_yaml_all_testfiles
    else: 
        raise ValueError("specify if \'quicktest\' or \'all_tests\' on command line")

    testfieldownloader = alphaquant.benchm.testfile_handling.TestFileDownloader(test_folder=test_folder, links_yaml=links_yaml)
    testfieldownloader.download_missing_files()