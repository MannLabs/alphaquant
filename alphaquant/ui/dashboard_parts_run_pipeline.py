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
import alphaquant.run_pipeline as diffmgr
import alphaquant.ui.dashboad_parts_plots_basic as dashboad_parts_plots_basic

import alphabase.quantification.quant_reader.config_dict_loader as config_dict_loader
config_dict_loader.INTABLE_CONFIG = os.path.join(pathlib.Path(__file__).parent.absolute(), "../config/quant_reader_config_for_gui.yaml")

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
	def __init__(self, state, **params):
		super().__init__(**params)
		self.state = state
		self._setup_matplotlib()
		self._setup_logger()
		self._make_widgets()
		self.layout = None

	def _setup_matplotlib(self):
		"""Configure matplotlib to use a non-GUI backend and turn off interactive mode."""
		import matplotlib
		matplotlib.use('agg')
		import matplotlib.pyplot as plt
		plt.ioff()

	def _setup_logger(self):
		"""Configure logging to capture output in the console widget."""
		import logging
		import io

		# Create a StringIO object to capture log output
		self.log_stream = io.StringIO()

		# Create a handler that writes to our StringIO object
		self.stream_handler = logging.StreamHandler(self.log_stream)
		self.stream_handler.setLevel(logging.INFO)

		# Create a formatter and set it for the handler
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		self.stream_handler.setFormatter(formatter)

		# Get the root logger and add our handler
		self.logger = logging.getLogger()
		self.logger.addHandler(self.stream_handler)
		self.logger.setLevel(logging.INFO)

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

		# Create a Row with the Select widget and its description
		self.sample_mapping_select = pn.widgets.Select(
			name='Sample Mapping Mode',
			options=[
				'Upload sample to condition file',
				'Generate new sample to condition map'
			],
			value='Upload sample to condition file',
			width=300,
			description='Choose whether to upload an existing sample-to-condition mapping file or create a new one'
		)

		self.sample_mapping_mode_container = pn.Row(
			self.sample_mapping_select,
			align='start'
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
			name='Organism',
			options=['human', 'mouse'],
			value='human',
			width=300,
			description='Select the organism your samples come from'
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

		# Update threshold settings widgets with descriptions
		self.minrep_either = pn.widgets.IntInput(
			name='Min replicates (either condition):',
			value=2,
			start=0,
			width=300,
			description='Minimum number of valid values required in at least one of the conditions'
		)

		self.minrep_both = pn.widgets.IntInput(
			name='Min replicates (both conditions):',
			value=2,
			start=1,
			width=300,
			description='Minimum number of valid values required in both conditions'
		)

		self.minrep_c1 = pn.widgets.IntInput(
			name='Min replicates (condition 1):',
			value=2,
			start=0,
			width=300,
			description='Minimum number of valid values required in condition 1',
			visible=False
		)

		self.minrep_c2 = pn.widgets.IntInput(
			name='Min replicates (condition 2):',
			value=2,
			start=0,
			width=300,
			description='Minimum number of valid values required in condition 2',
			visible=False
		)

		self.min_num_ions = pn.widgets.IntInput(
			name='Min number of ions per peptide:',
			value=1,
			start=1,
			width=300,
			description='Minimum number of ions required for each peptide to be included in the analysis'
		)

		self.minpep = pn.widgets.IntInput(
			name='Min peptides per protein:',
			value=1,
			start=1,
			width=300,
			description='Minimum number of peptides required for each protein to be included in the analysis'
		)

		self.cluster_threshold_pval = pn.widgets.FloatInput(
			name='Clustering p-value threshold:',
			value=0.001,
			start=0,
			end=1,
			width=300,
			description='P-value threshold used for clustering analysis'
		)

		self.volcano_fdr = pn.widgets.FloatInput(
			name='Volcano plot FDR:',
			value=0.05,
			start=0,
			end=1,
			width=300,
			description='False Discovery Rate threshold for the volcano plot'
		)

		self.volcano_fcthresh = pn.widgets.FloatInput(
			name='Volcano plot fold change threshold:',
			value=0.5,
			start=0,
			width=300,
			description='Fold change threshold for highlighting significant changes in the volcano plot'
		)

		self.condition_comparison_header = pn.pane.Markdown(
		"### Available Condition Comparisons",
		visible=True
		)

		self.condition_comparison_instructions = pn.pane.Markdown(
			"Select the condition pairs you want to analyze:",
			visible=True
		)

		# Replace the analysis_type Select widget
		self.analysis_type = pn.widgets.Select(
			name='Select Condition Analysis Type',
			options=['Pairwise Comparison', 'Median Condition Analysis'],
			value='Select an analysis',
			description='Choose between comparing pairs of conditions or comparing each condition against a median reference'
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
		self.run_pipeline_error = pn.pane.Alert(
			alert_type="danger",
			visible=False,
			margin=(5, 5, 5, 5),
		)

		# Add console output widget
		self.console_output = pn.widgets.TextAreaInput(
			placeholder='Pipeline output will appear here...',
			height=800,
			width=400,
			disabled=True
		)

		# Replace the visualize_data_button with a progress panel
		self.condition_progress_panel = pn.Column(
			pn.pane.Markdown("### Analysis Progress", margin=(0,0,10,0)),
			pn.pane.Markdown(
				"As soon as a condition pair is finished, you can inspect the results in the plots tabs.",
				margin=(0,0,10,0)
			),
			sizing_mode='stretch_width'
		)

		# Dictionary to store progress indicators for each condition pair
		self.condition_progress = {}

		# Watchers
		self.sample_mapping_select.param.watch(self._toggle_sample_mapping_mode, 'value')
		self.path_analysis_file.param.watch(
			self._activate_after_analysis_file_upload, 'value'
		)
		self.samplemap_fileupload.param.watch(self._update_samplemap_table, 'value')
		self.samplemap_table.param.watch(self._add_conditions_for_assignment, 'value')
		self.minrep_either.param.watch(self._update_minrep_both, 'value')
		self.run_pipeline_button.param.watch(self._run_pipeline, 'clicks')
		self.analysis_type.param.watch(self._toggle_analysis_type, 'value')
		self.filtering_options.param.watch(self._toggle_filtering_options, 'value')
		self.path_output_folder.param.watch(self._update_results_dir, 'value')
		self.path_analysis_file.param.watch(self._update_analysis_file, 'value')
		self.samplemap_fileupload.param.watch(self._update_samplemap, 'value')


	def create(self):
		"""
		Build and return the main layout for the pipeline widget.
		"""
		# Advanced Configuration Card
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
			self.sample_mapping_mode_container,
			pn.Row(
				self.loading_samples_indicator,
				self.loading_samples_message
			),
			self.samplemap_fileupload,
			self.samplemap_table
		)

		# Create condition comparison layout
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

		# Main column without scroll
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
				self.run_pipeline_progress,
				sizing_mode='stretch_width'
			),
			self.run_pipeline_error,
			self.condition_progress_panel,  # Add progress panel here
			sizing_mode='stretch_width'
		)

		# Console output column wrapped in a Row for padding
		console_col = pn.Row(
			pn.Column(
				self.console_output,
				width=400
			),
			width=425,
			align='start'
		)

		# Main layout with frame
		self.layout = pn.Column(
			pn.Row(
				main_col,
				console_col,
				sizing_mode='stretch_width'
			),
			css_classes=['custom-frame'],  # For custom styling if needed
			margin=(20, 20, 20, 20),  # This handles both outer and inner spacing
			sizing_mode='stretch_width',
			styles={
				'background': '#f8f9fa',    # Light gray background
				'border': '1px solid #dee2e6',  # Light gray border
				'border-radius': '5px',         # Rounded corners
				'box-shadow': '0 1px 3px rgba(0,0,0,0.12)'  # Subtle shadow
			}
		)

		return self.layout

	def _run_pipeline(self, *events):
		"""Run the alphaquant pipeline when the button is clicked."""
		if self.analysis_type.value == 'Select an analysis':
			self.run_pipeline_error.object = "Please select an analysis type before running the pipeline."
			self.run_pipeline_error.visible = True
			return

		self.run_pipeline_progress.active = True
		self.run_pipeline_error.visible = False
		self.console_output.value = "Starting pipeline...\n"

		# Update progress panel with selected condition pairs
		self._update_progress_panel(self.assign_cond_pairs.value or [])

		# Initial check for existing results and start monitoring
		self._check_existing_results()
		if not hasattr(self, 'results_monitor'):
			self.results_monitor = pn.state.add_periodic_callback(
				self._check_existing_results,
				period=1000  # Check every second
			)

		try:
			# Get condition combinations
			if self.assign_cond_pairs.value:
				cond_combinations = [
					tuple(pair.split('_VS_'))
					for pair in self.assign_cond_pairs.value
				]
			else:
				cond_combinations = [
					tuple(pair.split('_VS_'))
					for pair in self.assign_cond_pairs.options
				]

			# Collect all configuration parameters
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

			# Update console output
			self.console_output.value += f"Running pipeline with parameters:\n{pipeline_params}\n"

			# Run the pipeline
			diffmgr.run_pipeline(**pipeline_params)

			self.logger.info("Pipeline completed successfully!")

		except Exception as e:
			error_message = f"Error running pipeline: {e}"
			self.run_pipeline_error.object = error_message
			self.run_pipeline_error.visible = True
			self.logger.error(error_message)

		finally:
			# Stop the monitoring
			if hasattr(self, 'results_monitor'):
				self.results_monitor.stop()
				delattr(self, 'results_monitor')

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


	def _activate_after_analysis_file_upload(self, event):
		"""Handle analysis file upload."""
		if event.new:
			self._set_default_output_folder()
			self.path_output_folder.disabled = False
			self.run_pipeline_button.disabled = False

	def _set_default_output_folder(self):
		"""Set default output folder and start monitoring."""
		if self.path_analysis_file.value:
			base_path = os.path.dirname(self.path_analysis_file.value)
			output_path = os.path.join(base_path, 'results')
			self.path_output_folder.value = output_path
			self.state.results_dir = output_path

			# Make sure we have a progress panel
			if not hasattr(self, 'condition_progress'):
				self._update_progress_panel(self.assign_cond_pairs.value or [])

			# Start monitoring if not already running
			if not hasattr(self, 'results_monitor'):
				self.results_monitor = pn.state.add_periodic_callback(
					self._check_existing_results,
					period=1000
				)

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
		When a sample map file is uploaded, parse it into a DataFrame and update the state.
		"""
		if not self.samplemap_fileupload.value:
			return

		file_ext = os.path.splitext(self.samplemap_fileupload.filename)[-1].lower()
		sep = ',' if file_ext == '.csv' else '\t'

		try:
			# Parse the uploaded file into DataFrame
			df = pd.read_csv(
				StringIO(self.samplemap_fileupload.value.decode('utf-8')),
				sep=sep,
				dtype=str
			)

			# Update the table widget
			self.samplemap_table.value = df
			self.samplemap_table.visible = True

			# Update the state with the DataFrame
			self.state.samplemap_df = df
			self.state.notify_subscribers('samplemap_df')

		except Exception as e:
			self.run_pipeline_error.object = f"Error reading sample map: {e}"
			self.run_pipeline_error.visible = True

	def _add_conditions_for_assignment(self, *events):
		"""Update conditions and immediately check for existing results."""
		if self.samplemap_table.value is None:
			return

		df = self.samplemap_table.value
		if 'condition' in df.columns:
			unique_condit = df['condition'].dropna().unique()
			comb_condit = [
				'_VS_'.join(comb)
				for comb in itertools.permutations(unique_condit, 2)
			]
			self.assign_cond_pairs.options = comb_condit

			# Update progress panel and check for existing results
			self._update_progress_panel(comb_condit)
			self._check_existing_results()

	def _update_progress_panel(self, condition_pairs):
		"""Simplified progress panel creation."""
		selected_pairs = self.assign_cond_pairs.value or []

		# Clear existing progress indicators
		self.condition_progress.clear()
		progress_items = []

		# Overall status
		self.overall_status = pn.pane.Markdown(
			"### Analysis Progress",
			margin=(0, 0, 10, 0)
		)
		progress_items.append(self.overall_status)

		for pair in selected_pairs:
			# Simple status row with all indicators
			spinner = pn.indicators.LoadingSpinner(value=False, color='primary', width=20, height=20)
			pending = pn.pane.Markdown('⏳', styles={'font-size': '20px'})
			running = pn.pane.Markdown('🔄', styles={'font-size': '20px'})
			complete = pn.pane.Markdown('✅', styles={'color': 'green', 'font-size': '20px'})

			# Set initial visibility
			spinner.visible = False
			pending.visible = True
			running.visible = False
			complete.visible = False

			status_row = pn.Row(
				spinner,
				pending,
				running,
				complete,
				pn.pane.Markdown(f"**{pair}**"),
				margin=(5, 5, 5, 10)
			)

			self.condition_progress[pair] = {
				'row': status_row,
				'spinner': spinner,
				'pending': pending,
				'running': running,
				'complete': complete
			}
			progress_items.append(status_row)

		# Create the card
		progress_card = pn.Card(
			pn.Column(*progress_items, margin=10),
			margin=5
		)

		# Update panel
		self.condition_progress_panel.clear()
		self.condition_progress_panel.append(progress_card)

	def _check_existing_results(self):
		"""Simple check for existing result files."""
		output_dir = self.path_output_folder.value
		if not output_dir or not os.path.exists(output_dir):
			return

		print(f"\nChecking results in: {output_dir}")

		for pair in self.condition_progress:
			cond1, cond2 = pair.split('_VS_')
			result_file = os.path.join(output_dir, f"{cond1}_VS_{cond2}.results.tsv")

			print(f"Checking file: {result_file}")
			file_exists = os.path.exists(result_file)
			print(f"File exists: {file_exists}")

			if file_exists:
				print(f"Setting indicators for {pair}")
				try:
					self.condition_progress[pair]['pending'].visible = False
					self.condition_progress[pair]['running'].visible = False
					self.condition_progress[pair]['spinner'].visible = False
					self.condition_progress[pair]['complete'].visible = True
					print(f"Successfully marked {pair} as complete")
				except Exception as e:
					print(f"Error updating indicators for {pair}: {str(e)}")
			else:
				self.condition_progress[pair]['pending'].visible = True
				self.condition_progress[pair]['running'].visible = False
				self.condition_progress[pair]['spinner'].visible = False
				self.condition_progress[pair]['complete'].visible = False

	def _update_minrep_both(self, *events):
		"""Set minrep_both to 0 when minrep_either is changed."""
		self.minrep_both.value = 0

	def _update_results_dir(self, event):
		"""Update central state with new results directory and check for existing files."""
		value = event
		if isinstance(event, param.Event):
			value = event.new

		if value:
			print(f"\nOutput directory changed to: {value}")
			self.state.results_dir = str(value)

			# Create progress panel if it doesn't exist
			if not hasattr(self, 'condition_progress'):
				self._update_progress_panel(self.assign_cond_pairs.value or [])

			# Immediate check for existing files
			self._check_existing_results()

			# Start monitoring if not already running
			if not hasattr(self, 'results_monitor'):
				self.results_monitor = pn.state.add_periodic_callback(
					self._check_existing_results,
					period=1000
				)

	def _update_analysis_file(self, event):
		"""Update central state with new analysis file."""
		self.state.analysis_file = event.new
		self.state.notify_subscribers('analysis_file')

	def _update_samplemap(self, event):
		"""Update central state with new sample map file."""
		if event.new:
			try:
				# Parse the uploaded file into DataFrame
				file_ext = os.path.splitext(self.samplemap_fileupload.filename)[-1].lower()
				sep = ',' if file_ext == '.csv' else '\t'

				df = pd.read_csv(
					StringIO(event.new.decode('utf-8')),
					sep=sep,
					dtype=str
				)

				# Update the state with the DataFrame
				self.state.samplemap_df = df
				self.state.notify_subscribers('samplemap_df')

			except Exception as e:
				self.run_pipeline_error.object = f"Error reading sample map: {str(e)}"
				self.run_pipeline_error.visible = True

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

	def _update_console(self):
		"""Update the console output widget with new log messages."""
		# Get any new log messages
		self.log_stream.seek(0)
		new_logs = self.log_stream.read()

		# Update the console widget
		if new_logs:
			current_text = self.console_output.value or ""
			self.console_output.value = current_text + new_logs

			# Clear the StringIO buffer
			self.log_stream.truncate(0)
			self.log_stream.seek(0)

	def __del__(self):
		"""Clean up logging handler when the widget is destroyed."""
		if hasattr(self, 'stream_handler'):
			self.logger.removeHandler(self.stream_handler)
			self.stream_handler.close()

class Tabs(param.Parameterized):
	"""
	This class creates a single pn.Tabs layout containing:
	  1. Pipeline
	  2. Single Comparison
	  3. Plotting
	"""
	pipeline = param.ClassSelector(class_=RunPipeline)
	main_tabs = param.ClassSelector(class_=pn.Tabs, allow_None=True)

	def __init__(self, pipeline, **params):
		super().__init__(pipeline=pipeline, **params)
		self._build_initial_tabs()
		# Watch for changes in the pipeline's output folder
		self.pipeline.path_output_folder.param.watch(self._sync_results_dir, 'value')
		self.pipeline.param.watch(
			self._update_tabs,
			'update_event',
			onlychanged=True
		)

	def _sync_results_dir(self, event):
		"""Sync the results directory between Pipeline and Plotting tabs."""
		if event.new:  # Only update if there's a value
			try:
				# Update the Plotting tab's results directory
				plotting_tab = self.main_tabs[1][1]  # Access the plotting tab content
				if isinstance(plotting_tab, dashboad_parts_plots_basic.PlottingTab):
					plotting_tab.results_dir_input.value = event.new
			except Exception as e:
				print(f"Error syncing results directory: {str(e)}")

	def _build_initial_tabs(self):
		"""Create initial empty tabs."""
		self.main_tabs = pn.Tabs(
			('Single Comparison', pn.pane.Markdown(
				"## No data loaded\nPlease load data in the Pipeline tab first."
			)),
			('Plotting', dashboad_parts_plots_basic.PlottingTab().panel()),
			tabs_location='above',
			sizing_mode='stretch_width',
			margin=(10, 10, 10, 10)
		)

	def _update_tabs(self, event=None):
		"""Update tabs with visualization when data is available."""
		try:
			if (self.pipeline.path_output_folder.value and
				self.pipeline.samplemap_table.value is not None):

				# Update Plotting tab
				plotting_tab = dashboad_parts_plots_basic.PlottingTab(
					results_dir=self.pipeline.path_output_folder.value
				)
				self.main_tabs[1] = ('Plotting', plotting_tab.panel())

		except Exception as e:
			error_msg = f"Error updating visualization tabs: {str(e)}"
			self.main_tabs[0] = ('Single Comparison', pn.pane.Markdown(
				f"### Visualization Error\n\n{error_msg}"
			))
			self.main_tabs[1] = ('Plotting', pn.pane.Markdown(
				f"### Visualization Error\n\n{error_msg}"
			))


def build_dashboard():
	"""Build the overall dashboard layout."""
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

	# Create pipeline instance
	pipeline = RunPipeline(state=None)
	pipeline_layout = pipeline.create()

	# Create tab manager with pipeline tab
	tab_manager = Tabs(pipeline)
	tabs = tab_manager.create()

	# Create tabs with Pipeline as the first tab
	all_tabs = pn.Tabs(
		('Pipeline', pipeline_layout),
		('Single Comparison', tabs[0][1]),
		('Plotting', tabs[1][1]),
		dynamic=True,
		tabs_location='above',
		sizing_mode='stretch_width'
	)

	# Main layout
	main_layout = pn.Column(
		header.create(),
		pn.layout.Divider(),
		main_text.create(),
		all_tabs,
		sizing_mode='stretch_width'
	)

	template = pn.template.FastListTemplate(
		title="AlphaQuant Analysis",
		sidebar=[],
		main=[main_layout],
		theme='dark',
		main_max_width="1200px",
		main_layout="width"
	)
	return template

# If run as a script, you can do:
# if __name__ == "__main__":
#     dash = build_dashboard()
#     dash.servable()
