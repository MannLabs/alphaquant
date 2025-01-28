import panel as pn

class PlottingTab():
	def __init__(self, results_dir_phospho, results_dir_proteome, cond_pairs):
		self.results_dir_phospho = results_dir_phospho
		self.results_dir_proteome = results_dir_proteome
		self.cond_pairs = cond_pairs

	def panel(self):
		return pn.pane.Markdown(f"## Plotting Tab")
