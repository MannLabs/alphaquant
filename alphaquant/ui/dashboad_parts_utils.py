import os

def get_samplemap_file_path(results_dir):
	"""Get the appropriate samplemap file path based on the results directory."""
	if os.path.exists(os.path.join(results_dir, 'samplemap_w_median.tsv')):
		return os.path.join(results_dir, 'samplemap_w_median.tsv')
	else:
		return os.path.join(results_dir, 'samplemap.tsv')
