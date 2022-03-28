import alphaquant.testfile_handling
test_folder = "."
links_yaml = "../alphaquant/configs/download_links_for_testfiles.yaml"


testfieldownloader = alphaquant.testfile_handling.TestFileDownloader(test_folder=test_folder, links_yaml=links_yaml)
testfieldownloader.download_missing_files()