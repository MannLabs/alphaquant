import os
import re
from io import StringIO
import itertools
import pathlib

import param
import panel as pn
import pandas as pd
import matplotlib
matplotlib.use('agg')

# alphaquant imports
import alphaquant.diffquant.diffutils as aqdiffutils
import alphaquant.run_pipeline as diffmgr
import alphaquant.ui.gui_textfields as gui_textfields

import alphabase.quantification.quant_reader.config_dict_loader as config_dict_loader
config_dict_loader.INTABLE_CONFIG = os.path.join(pathlib.Path(__file__).parent.absolute(), "../config/quant_reader_config.yaml")

# If using Plotly in Panel
pn.extension('plotly')


class BaseWidget(param.Parameterized):
	"""
	A base class to hold a common event trigger
	for Panel/Param watchers in child classes.
	"""
	name = param.String()
	update_event = param.Integer(default=0, doc="Triggers re-computation")

	def __init__(self, name='', **params):
		super().__init__(name=name, **params)

	def trigger_dependency(self):
		"""Increment update_event to force watchers to re-run."""
		self.update_event += 1


class HeaderWidget(param.Parameterized):
	"""
	Header layout with project logos and links.
	"""
	title = param.String()
	img_folder_path = param.String()
	github_url = param.String()

	def __init__(self, title, img_folder_path, github_url, **params):
		super().__init__(
			title=title,
			img_folder_path=img_folder_path,
			github_url=github_url,
			**params
		)
		self._create_panes()

	def _create_panes(self):
		"""
		Initialize all Pane objects used in the header.
		"""
		self.header_title = pn.pane.Markdown(
			f'# {self.title}',
			sizing_mode='stretch_width'
		)
		self.biochem_logo_path = os.path.join(self.img_folder_path, "mpi_logo.png")
		self.mpi_logo_path = os.path.join(self.img_folder_path, "max-planck-gesellschaft.jpg")
		self.github_logo_path = os.path.join(self.img_folder_path, "github.png")

		# Logos (png & jpg) with links
		self.mpi_biochem_logo = pn.pane.PNG(
			self.biochem_logo_path,
			link_url='https://www.biochem.mpg.de/mann',
			width=60, height=60
		)
		self.mpi_logo = pn.pane.JPG(
			self.mpi_logo_path,
			link_url='https://www.biochem.mpg.de/en',
			height=60, width=60
		)
		self.github_logo = pn.pane.PNG(
			self.github_logo_path,
			link_url=self.github_url,
			height=60
		)

	def create(self):
		"""
		Return a layout (Panel Row) of the header widgets.
		"""
		return pn.Row(
			self.mpi_biochem_logo,
			self.mpi_logo,
			self.header_title,
			self.github_logo,
			sizing_mode='stretch_width',
			margin=(5, 10, 5, 10)
		)


class MainWidget(param.Parameterized):
	"""
	Create a layout for tool description and an optional manual download.
	"""
	description = param.String()
	manual_path = param.String()

	def __init__(self, description, manual_path, **params):
		super().__init__(description=description, manual_path=manual_path, **params)
		self._create_widgets()

	def _create_widgets(self):
		"""
		Initialize your markdown description and file downloader.
		"""
		self.project_description = pn.pane.Markdown(
			self.description,
			margin=(10, 0, 10, 0),
			sizing_mode='stretch_width'
		)

	def create(self):
		"""
		Return a simple column layout with the main description (and manual download).
		"""
		layout = pn.Column(
			self.project_description,
			sizing_mode='stretch_width',
			margin=(5, 5, 5, 5)
		)
		return layout


class RunPipeline(BaseWidget):
	"""
	Widget to gather file inputs, define condition pairs, and run an analysis pipeline.
	Includes advanced configuration options.
	"""
	def __init__(self, **params):
		super().__init__(**params)
		self._setup_matplotlib()
		self._make_widgets()
		self.layout = None

	def _setup_matplotlib(self):
		"""Configure matplotlib to use a non-GUI backend and turn off interactive mode."""
		import matplotlib
		matplotlib.use('agg')
		import matplotlib.pyplot as plt
		plt.ioff()

	def _make_widgets(self):
		"""
		Create all Panel/Param widgets used in the layout.
		"""
		# Add a new loading indicator
		self.loading_samples_indicator = pn.indicators.LoadingSpinner(
			value=False,
			color='primary',
			visible=False
		)
		self.loading_samples_message = pn.pane.Markdown(
			"Loading sample names...",
			visible=False
		)

		# File paths with descriptions
		self.path_analysis_file = pn.widgets.TextInput(
			name='Analysis file:',
			placeholder='Path to MQ/Spectronaut/DIA-NN file',
			width=700,
			sizing_mode='fixed',
			description='Enter the full path to your input file from MaxQuant, Spectronaut, or DIA-NN'
		)
		self.path_output_folder = pn.widgets.TextInput(
			name='Output folder:',
			placeholder='Path to output folder',
			width=700,
			sizing_mode='fixed',
			description='Specify where you want the analysis results to be saved'
		)

		self.sample_mapping_mode = pn.widgets.Select(
			options=[
				'Upload sample to condition file',
				'Generate new sample to condition map'
			],
			value='Upload sample to condition file',
			width=300,
			description='Choose whether to upload an existing sample-to-condition mapping file or create a new one'
		)

		self.samplemap_fileupload = pn.widgets.FileInput(
			accept='.tsv,.csv,.txt',
			margin=(5, 5, 10, 20),
			visible=True  # Add this line explicitly
		)
		# In _make_widgets(), when initializing samplemap_table, add visible=False:
		self.samplemap_table = pn.widgets.Tabulator(
			layout='fit_data_fill',
			height=250,
			show_index=False,
			width=500,
			margin=(5, 5, 5, 5),
			visible=False  # Add this line
		)
		self.assign_cond_pairs = pn.widgets.CrossSelector(
			width=600,
			height=250,
			margin=(5, 5, 10, 20),
			name='Select condition pairs'
		)

		# Advanced configuration widgets with descriptions
		self.modification_type = pn.widgets.TextInput(
			name='Modification type:',
			placeholder='e.g., [Phospho (STY)] for Spectronaut',
			width=300,
			description='Specify the modification type exactly as it appears in your data (e.g., [Phospho (STY)] for Spectronaut)'
		)
		self.input_type = pn.widgets.TextInput(
			name='Input type:',
			placeholder='Type of quantitative information',
			width=300
		)
		self.organism = pn.widgets.Select(
			name='Organism:',
			options=['human', 'mouse'],
			value='human',  # Set default value
			width=300
		)
		self.filtering_options = pn.widgets.Select(
			name='Filtering Options:',
			options=[
				'min. valid values in condition1 OR condition2',
				'min. valid values in condition1 AND condition2',
				'set min. valid values per condition'
			],
			value='min. valid values in condition1 OR condition2',
			width=300,
			description='Choose how to filter your data based on the number of valid values in each condition'
		)

		self.minrep_either = pn.widgets.IntInput(
			value=2,
			start=0,
			width=300,
			visible=True  # Initially visible as it's the default option
		)

		self.minrep_both = pn.widgets.IntInput(
			value=2,
			start=1,
			width=300,
			visible=False
		)

		self.minrep_c1 = pn.widgets.IntInput(
			name='Min replicates (condition 1):',
			value=2,
			start=0,
			width=300,
			visible=False
		)

		self.minrep_c2 = pn.widgets.IntInput(
			name='Min replicates (condition 2):',
			value=2,
			start=0,
			width=300,
			visible=False
		)

		self.min_num_ions = pn.widgets.IntInput(
			name='Min number of ions per peptide:',
			value=1,
			start=1,
			width=300
		)
		self.minpep = pn.widgets.IntInput(
			name='Min peptides per protein:',
			value=1,
			start=1,
			width=300
		)
		self.cluster_threshold_pval = pn.widgets.FloatInput(
			name='Clustering p-value threshold:',
			value=0.001,
			start=0,
			end=1,
			width=300
		)
		self.volcano_fdr = pn.widgets.FloatInput(
			name='Volcano plot FDR:',
			value=0.05,
			start=0,
			end=1,
			width=300
		)
		self.volcano_fcthresh = pn.widgets.FloatInput(
			name='Volcano plot fold change threshold:',
			value=0.5,
			start=0,
			width=300
		)

		self.condition_comparison_header = pn.pane.Markdown(
		"### Available Condition Comparisons",
		visible=True
		)

		self.condition_comparison_instructions = pn.pane.Markdown(
			"Select the condition pairs you want to analyze:",
			visible=True
		)

		# Replace the medianref_analysis_switch with a dropdown menu
		self.analysis_type = pn.widgets.Select(
			name='Select Condition Analysis Type:',
			options=['Pairwise Comparison', 'Median Condition Analysis'],
			value='Select an analysis'
		)

		# A pane for showing the "comparing every condition..." message
		# which is hidden by default
		self.medianref_message = pn.pane.Markdown(
			"Every condition will be compared against the median reference",
			visible=False,  # start hidden
		)
		# Boolean switches with descriptions
		self.switches = {
			'use_ml': pn.widgets.Switch(
				name='Enable machine learning',
				value=True
			),
			'take_median_ion': pn.widgets.Switch(
				name='Use median-centered ions',
				value=True
			),
			'perform_ptm_mapping': pn.widgets.Switch(
				name='Enable PTM mapping',
				value=False
			),
			'perform_phospho_inference': pn.widgets.Switch(
				name='Enable phospho inference',
				value=False
			),
			'outlier_correction': pn.widgets.Switch(
				name='Enable outlier correction',
				value=True
			),
			'normalize': pn.widgets.Switch(
				name='Enable normalization',
				value=True
			),
			'use_iontree_if_possible': pn.widgets.Switch(
				name='Use ion tree when possible',
				value=True
			),
			'write_out_results_tree': pn.widgets.Switch(
				name='Write results tree',
				value=True
			),
			'use_multiprocessing': pn.widgets.Switch(
				name='Enable multiprocessing',
				value=False
			),
			'runtime_plots': pn.widgets.Switch(
				name='Generate runtime plots',
				value=True
			),
		}

		# If you want to keep the descriptions, you can add them separately as Markdown panes
		self.switch_descriptions = {
			'use_ml': pn.pane.Markdown('Use machine learning for improved data analysis'),
			'take_median_ion': pn.pane.Markdown('Center ion intensities around their median values'),
			'perform_ptm_mapping': pn.pane.Markdown('Map post-translational modifications to proteins'),
			'perform_phospho_inference': pn.pane.Markdown('Infer phosphorylation sites from the data'),
			'outlier_correction': pn.pane.Markdown('Automatically detect and correct outliers in the data'),
			'normalize': pn.pane.Markdown('Normalize data to account for technical variations'),
			'use_iontree_if_possible': pn.pane.Markdown('Use hierarchical ion structure when available'),
			'write_out_results_tree': pn.pane.Markdown('Save detailed results in a tree structure'),
			'use_multiprocessing': pn.pane.Markdown('Use multiple CPU cores to speed up processing (may use more memory)'),
			'runtime_plots': pn.pane.Markdown('Create plots during analysis to visualize the process'),
		}

		# Pipeline execution widgets with descriptions
		self.run_pipeline_button = pn.widgets.Button(
			name='Run pipeline',
			button_type='primary',
			height=35,
			width=170,
			margin=(10, 0, 0, 5),
			description='Start the analysis pipeline with current settings'
		)
		self.run_pipeline_progress = pn.indicators.Progress(
			active=False,
			bar_color='light',
			width=170,
			margin=(5, 0, 5, 5)
		)
		self.visualize_data_button = pn.widgets.Button(
			name='Visualize data',
			button_type='success',
			height=35,
			width=170,
			margin=(10, 0, 0, 5),
			description='View analysis results and generate visualizations'
		)
		self.run_pipeline_error = pn.pane.Alert(
			alert_type="danger",
			visible=False,
			margin=(5, 5, 5, 5),
		)

		# Watchers
		self.sample_mapping_mode.param.watch(self._toggle_sample_mapping_mode, 'value')
		self.path_analysis_file.param.watch(
			self._activate_after_analysis_file_upload, 'value'
		)
		self.samplemap_fileupload.param.watch(self._update_samplemap_table, 'value')
		self.samplemap_table.param.watch(self._add_conditions_for_assignment, 'value')
		self.minrep_either.param.watch(self._update_minrep_both, 'value')
		self.run_pipeline_button.param.watch(self._run_pipeline, 'clicks')
		self.visualize_data_button.param.watch(self._visualize_data, 'clicks')
		self.analysis_type.param.watch(self._toggle_analysis_type, 'value')
		self.filtering_options.param.watch(self._toggle_filtering_options, 'value')


	def create(self):
		"""
		Build and return the main layout for the pipeline widget.
		"""
		# 1) Instructions Card
		instructions_card = pn.Card(
			"### Instructions",
			gui_textfields.Descriptions.project_instruction,
			gui_textfields.Cards.spectronaut,
			gui_textfields.Cards.diann,
			gui_textfields.Cards.alphapept,
			gui_textfields.Cards.maxquant,
			title='Instructions',
			collapsed=True,
			margin=(5, 5, 5, 5),
			sizing_mode='fixed',
			width=400
		)


		# 2) Advanced Configuration Card
		advanced_settings_card = pn.Card(
			pn.Column(
				"### Threshold Settings",
				self.minrep_both,
				self.min_num_ions,
				self.minpep,
				self.cluster_threshold_pval,
				self.volcano_fdr,
				self.volcano_fcthresh,
				pn.layout.Divider(),
				"### Analysis Options",
				pn.Column(*[
					pn.Row(
						switch,
						self.switch_descriptions[key],
						align='center'
					) for key, switch in self.switches.items()
				]),
			),
			title='Advanced Configuration',
			collapsed=True,
			margin=(5, 5, 5, 5),
			sizing_mode='fixed',
			width=400
		)

		# Create samples and conditions layout
		samples_conditions_layout = pn.Column(
			self.sample_mapping_mode,
			pn.Row(
				self.loading_samples_indicator,
				self.loading_samples_message
			),
			self.samplemap_fileupload,
			self.samplemap_table
		)

		# Create condition comparison layout (without the header)
		condition_comparison_layout = pn.Column(
			self.condition_comparison_instructions,
			self.assign_cond_pairs,
			self.medianref_message,
		)

		# Create PTM settings card with fixed width
		ptm_settings_card = pn.Card(
			pn.Column(
				self.modification_type,
				self.organism,
			),
			title='PTM Settings',
			collapsed=True,
			margin=(5, 5, 5, 5),
			sizing_mode='fixed',
			width=400
		)

		# Main layout
		main_col = pn.Column(
			"### Input Files",
			self.path_analysis_file,
			self.path_output_folder,
			"### Samples and Conditions",
			samples_conditions_layout,
			self.analysis_type,
			condition_comparison_layout,
			"### Basic Settings",
			self.filtering_options,
			self.minrep_either,
			self.minrep_both,
			self.minrep_c1,
			self.minrep_c2,
			ptm_settings_card,
			advanced_settings_card,
			"### Pipeline Controls",
			pn.Row(
				self.run_pipeline_button,
				self.visualize_data_button,
				self.run_pipeline_progress,
				sizing_mode='stretch_width'
			),
			self.run_pipeline_error,
			sizing_mode='stretch_width'
		)

		# Show/hide components based on selected analysis type
		if self.analysis_type.value == 'Select an analysis':
			self.medianref_message.visible = False
			self.assign_cond_pairs.visible = False
			self.condition_comparison_header.visible = False
			self.condition_comparison_instructions.visible = False
		else:
			if self.analysis_type.value == 'Median Condition Analysis':
				self.medianref_message.visible = True
				self.assign_cond_pairs.visible = False
				self.condition_comparison_header.visible = False
				self.condition_comparison_instructions.visible = False
			else:
				self.medianref_message.visible = False
				self.assign_cond_pairs.visible = True
				self.condition_comparison_header.visible = True
				self.condition_comparison_instructions.visible = True

		# Main pipeline card
		main_pipeline_card = pn.Card(
			main_col,
			title='Run Pipeline',
			header_color='#333',
			header_background='#eaeaea',
			sizing_mode='stretch_width',
			margin=(10, 10, 10, 10)
		)

		# Final layout
		self.layout = pn.Column(
			instructions_card,
			main_pipeline_card
		)
		return self.layout

	def _run_pipeline(self, *events):
		"""
		Run the alphaquant pipeline when the button is clicked.
		"""
		if self.analysis_type.value == 'Select an analysis':
			self.run_pipeline_error.object = "Please select an analysis type before running the pipeline."
			self.run_pipeline_error.visible = True
			return


		self.run_pipeline_progress.active = True
		self.run_pipeline_error.visible = False

		# Get condition combinations from the CrossSelector
		if self.assign_cond_pairs.value:
			cond_combinations = [
				tuple(pair.split('_vs_'))
				for pair in self.assign_cond_pairs.value
			]
		else:
			cond_combinations = [
				tuple(pair.split('_vs_'))
				for pair in self.assign_cond_pairs.options
			]

		try:
			# Collect all configuration parameters
			print(self.samplemap_table.value)
			pipeline_params = {
				'input_file': self.path_analysis_file.value,
				'samplemap_df': self.samplemap_table.value,
				'results_dir': self.path_output_folder.value,
				'condpairs_list': cond_combinations,
				'modification_type': self.modification_type.value or None,
				'input_type_to_use': self.input_type.value or None,
				'organism': self.organism.value or None,
				'min_num_ions': self.min_num_ions.value,
				'minpep': self.minpep.value,
				'cluster_threshold_pval': self.cluster_threshold_pval.value,
				'volcano_fdr': self.volcano_fdr.value,
				'volcano_fcthresh': self.volcano_fcthresh.value,
			}
			# Add all switch values
			pipeline_params.update({
				key: switch.value for key, switch in self.switches.items()
			})

			if self.filtering_options.value == 'min. valid values in condition1 OR condition2':
				pipeline_params['minrep_either'] = self.minrep_either.value
			elif self.filtering_options.value == 'min. valid values in condition1 AND condition2':
				pipeline_params['minrep_both'] = self.minrep_both.value
			else:  # set min. valid values per condition
				pipeline_params['minrep_c1'] = self.minrep_c1.value
				pipeline_params['minrep_c2'] = self.minrep_c2.value


			# Run the pipeline
			diffmgr.run_pipeline(**pipeline_params)

		except Exception as e:
			self.run_pipeline_error.object = f"Error running pipeline: {e}"
			self.run_pipeline_error.visible = True

		self.trigger_dependency()
		self.run_pipeline_progress.active = False

		# Show/hide components based on selected analysis type
		if self.analysis_type.value == 'Median Condition Analysis':
			# Show components related to median condition analysis
			self.medianref_message.visible = True
			self.assign_cond_pairs.visible = False
			self.condition_comparison_header.visible = False
			self.condition_comparison_instructions.visible = False
		else:
			# Show components related to pairwise comparison
			self.medianref_message.visible = False
			self.assign_cond_pairs.visible = True
			self.condition_comparison_header.visible = True
			self.condition_comparison_instructions.visible = True

	def _toggle_sample_mapping_mode(self, event):
		"""Toggle visibility of sample mapping components based on selected mode."""
		if event.new == 'Upload sample to condition file':
			self.samplemap_fileupload.visible = True
			# Only show table if there's data in it
			#self.samplemap_table.visible = hasattr(self, 'data') and not self.data.empty
		else:
			self.samplemap_fileupload.visible = False
			self.samplemap_table.visible = True
			self._import_sample_names()
			self._init_samplemap_df_template()


	def _activate_after_analysis_file_upload(self, *events):
		"""
		When a new analysis file is entered, set default output folder and import data.
		"""
		self._set_default_output_folder()



	def _set_default_output_folder(self):
		if (not self.path_output_folder.value) and self.path_analysis_file.value:
			base_path = os.path.dirname(self.path_analysis_file.value)
			self.path_output_folder.value = os.path.join(base_path, 'results')

	def _import_sample_names(self):
		if self.path_analysis_file.value:
			try:
				# Show loading indicator
				self.loading_samples_indicator.visible = True
				self.loading_samples_message.visible = True

				input_file = self.path_analysis_file.value
				_, config_dict, sep = config_dict_loader.get_input_type_and_config_dict(input_file)
				sample_column = config_dict["sample_ID"]
				sample_names = set()
				for chunk in pd.read_csv(input_file, sep=sep, usecols=[sample_column], chunksize=400000):
					sample_names.update(chunk[sample_column].unique())
				self.sample_names = sample_names

			except Exception as e:
				self.run_pipeline_error.object = f"Error importing data: {e}"
				self.run_pipeline_error.visible = True
			finally:
				# Hide loading indicator
				self.loading_samples_indicator.visible = False
				self.loading_samples_message.visible = False

	def _init_samplemap_df_template(self):
		if hasattr(self, 'sample_names'):
			sample_names = self.sample_names
			sorted_names = self.natural_sort(sample_names)
			self.samplemap_table.value = pd.DataFrame({
				'sample': sorted_names,
				'condition': [''] * len(sorted_names)
			})

	def _update_samplemap_table(self, *events):
		"""
		When a sample map file is uploaded, parse it into the Tabulator widget and show the table.
		"""
		if not self.samplemap_fileupload.value:
			return
		file_ext = os.path.splitext(self.samplemap_fileupload.filename)[-1].lower()
		sep = ',' if file_ext == '.csv' else '\t'
		try:
			df = pd.read_csv(
				StringIO(self.samplemap_fileupload.value.decode('utf-8')),
				sep=sep,
				dtype=str
			)
			self.samplemap_table.value = df
			# Show the table after successful upload
			self.samplemap_table.visible = True
		except Exception as e:
			self.run_pipeline_error.object = f"Error reading sample map: {e}"
			self.run_pipeline_error.visible = True

	def _add_conditions_for_assignment(self, *events):
		"""
		Whenever the samplemap table is updated, auto-generate condition pairs.
		"""
		if self.samplemap_table.value is None:
			return
		df = self.samplemap_table.value
		if 'condition' in df.columns:
			unique_condit = df['condition'].dropna().unique()
			comb_condit = [
				'_vs_'.join(comb)
				for comb in itertools.permutations(unique_condit, 2)
			]
			self.assign_cond_pairs.options = comb_condit

	def _update_minrep_both(self, *events):
		"""Set minrep_both to 0 when minrep_either is changed."""
		self.minrep_both.value = 0

	def _visualize_data(self, *events):
		"""
		Trigger an update event for any dependent tabs/components to load results.
		"""
		self.run_pipeline_error.visible = False
		self.trigger_dependency()

	def natural_sort(self, l):
		"""
		Sort a list in a way that numbers are sorted numerically rather than alphabetically.
		"""
		convert = lambda text: int(text) if text.isdigit() else text.lower()
		alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
		return sorted(l, key=alphanum_key)

	def _toggle_analysis_type(self, event):
		"""Show/hide CrossSelector and message based on analysis type."""
		if event.new == 'Median Condition Analysis':
			# Hide CrossSelector and its explanatory text
			self.assign_cond_pairs.visible = False
			self.condition_comparison_header.visible = False
			self.condition_comparison_instructions.visible = False
			# Show median reference message
			self.medianref_message.visible = True
		else:
			self.assign_cond_pairs.visible = True
			self.condition_comparison_header.visible = True
			self.condition_comparison_instructions.visible = True
			self.medianref_message.visible = False

	def _toggle_filtering_options(self, event):
		"""Toggle visibility of replicate input fields based on filtering option."""
		# Hide all first
		self.minrep_either.visible = False
		self.minrep_both.visible = False
		self.minrep_c1.visible = False
		self.minrep_c2.visible = False

		# Show relevant widgets based on selection
		if event.new == 'min. valid values in condition1 OR condition2':
			self.minrep_either.visible = True
		elif event.new == 'min. valid values in condition1 AND condition2':
			self.minrep_both.visible = True
		else:  # set min. valid values per condition
			self.minrep_c1.visible = True
			self.minrep_c2.visible = True

class Tabs(param.Parameterized):
	"""
	Manages switching from the pipeline stage to one or more analysis tabs.
	"""
	pipeline = param.ClassSelector(class_=RunPipeline)
	main_tabs = param.ClassSelector(class_=pn.Tabs, allow_None=True)
	current_layout = param.Parameter()

	def __init__(self, pipeline, **params):
		super().__init__(pipeline=pipeline, **params)
		self.main_tabs = None
		self.current_layout = None
		# Watch pipeline.update_event to create the tabs after pipeline has run
		self.pipeline.param.watch(self._create_tabs_if_possible, 'update_event')

	def _create_tabs_if_possible(self, *events):
		"""
		If the user has clicked 'visualize data', show the analysis tabs.
		"""
		if (self.pipeline.path_output_folder.value
				and self.pipeline.visualize_data_button.clicks > 0):
			try:
				self._build_tabs()
				# Update the current layout
				self.current_layout = self.main_tabs
			except Exception as e:
				print(f"Error creating visualization tabs: {str(e)}")
				self.current_layout = pn.pane.Markdown(
					f"Error creating visualization: {str(e)}"
				)

	def _build_tabs(self):
		"""
		Build a Panel Tabs layout with visualization components.
		"""
		if self.main_tabs is None:
			self.main_tabs = pn.Tabs(
				tabs_location='above',
				sizing_mode='stretch_width',
				margin=(10, 10, 10, 10)
			)

			# Initialize visualization components
			try:
				# Single Comparison tab
				single_comp = SingleComparison(
					self.pipeline.path_output_folder.value,
					self.pipeline.samplemap_table.value
				)
				self.main_tabs.append(
					('Single Comparison', single_comp.layout)
				)
				# Could add more tabs here for multiple comparisons, etc.
			except Exception as e:
				print(f"Error initializing comparison tab: {str(e)}")
				self.main_tabs.append(
					('Error', pn.pane.Markdown(f"Error creating visualization: {str(e)}"))
				)

	def create(self):
		"""
		Return the current view - either tabs or a message.
		"""
		if self.current_layout is not None:
			return self.current_layout
		else:
			return pn.pane.Markdown(
				"## No results yet.\nPlease run the pipeline and click **Visualize data**.",
				sizing_mode='stretch_width'
			)

# ----------------
# BUILD DASHBOARD
# ----------------
def build_dashboard():
	"""
	Example function to build the overall dashboard layout in a FastListTemplate.
	"""
	header = HeaderWidget(
		title="AlphaQuant Dashboard",
		img_folder_path="./assets",
		github_url="https://github.com/<my_repo>"
	)
	main_text = MainWidget(
		description=(
			"Welcome to our analysis dashboard. "
			"Please load your data and run the pipeline."
		),
		manual_path="path/to/manual.pdf"
	)
	pipeline = RunPipeline()
	tab_manager = Tabs(pipeline)

	template = pn.template.FastListTemplate(
		title="AlphaQuant Analysis",
		sidebar=[   # If you want a sidebar, you can put items here
			# "## Sidebar Title",
			# pn.widgets.Select(options=['Item1','Item2']),
		],
		main=[
			header.create(),
			pn.layout.Divider(),
			main_text.create(),
			pipeline.create(),
			tab_manager.create()
		],
		theme='dark',          # or 'default'
		main_max_width="1200px"
	)
	return template

# If run as a script, you can do:
# if __name__ == "__main__":
#     dash = build_dashboard()
#     dash.servable()
