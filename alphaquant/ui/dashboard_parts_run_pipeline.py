import os
import re
from io import StringIO
import itertools
import pathlib
import time

import param
import panel as pn
import pandas as pd
import matplotlib
matplotlib.use('agg')

# alphaquant imports
import alphaquant.run_pipeline as diffmgr
import alphaquant.config.variables as aq_variables
import alphaquant.ui.dashboad_parts_plots_basic as dashboad_parts_plots_basic
import alphaquant.ui.dashboard_parts_plots_proteoforms as dashboad_parts_plots_proteoforms
import alphaquant.ui.gui as gui
import alphaquant.ui.gui_textfields as gui_textfields

import alphabase.quantification.quant_reader.config_dict_loader as config_dict_loader
config_dict_loader.INTABLE_CONFIG = os.path.join(pathlib.Path(__file__).parent.absolute(), "../config/quant_reader_config.yaml")
# If using Plotly in Panel
pn.extension('plotly')

MEDIAN_CONDITION_ANALYSIS = "Median Condition Analysis"
PAIRWISE_COMPARISON = "Pairwise Comparison"

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
		self._progress_monitor = None
		self._path_output_folder = None
		self._condition_progress = {}
		self._overall_status = None
		self._log_stream = None
		self._stream_handler = None
		self._logger = None

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
		self._log_stream = io.StringIO()

		# Create a handler that writes to our StringIO object
		self._stream_handler = logging.StreamHandler(self._log_stream)
		self._stream_handler.setLevel(logging.INFO)

		# Create a formatter and set it for the handler
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		self._stream_handler.setFormatter(formatter)

		# Get the root logger and add our handler
		self._logger = logging.getLogger()
		self._logger.addHandler(self._stream_handler)
		self._logger.setLevel(logging.INFO)

	def _make_widgets(self):
		"""
		Create all Panel/Param widgets used in the layout.
		"""
		self.loading_samples_indicator = pn.indicators.LoadingSpinner(
			value=False,
			color='primary',
			visible=False
		)
		self.loading_samples_message = pn.pane.Markdown(
			"Loading sample names... In case you have a large file, this can take some minutes.",
			visible=False
		)

		self.template_success_message = pn.pane.Alert(
			"Template has been generated. Please fill out the condition column in the table below. The template has also been saved to your output folder if you prefer to edit it with Excel or other applications.",
			alert_type="success",
			visible=False
		)

		self.path_analysis_file = pn.widgets.TextInput(
			name='Analysis file:',
			placeholder='Path to AlphaDIA/DIA-NN/Spectronaut/Fragpipe etc. file',
			width=700,
			sizing_mode='fixed',
			description=gui_textfields.Descriptions.tooltips['file_input']
		)
		self.path_output_folder = pn.widgets.TextInput(
			name='Output folder:',
			placeholder='Path to output folder',
			width=700,
			sizing_mode='fixed',
			description='Specify where you want the analysis results to be saved'
		)

		self.sample_mapping_select = pn.widgets.Select(
			name='Sample Mapping Mode',
			options=['Generate new sample to condition map', 'Upload sample to condition file'],
			value='Generate new sample to condition map',
			width=300,
			description=gui_textfields.Descriptions.tooltips['sample_mapping']
		)

		self.generate_samplemap_button = pn.widgets.Button(
			name='Generate Samplemap Template',
			button_type='primary',
			width=100,
			disabled=True,  # Initially disabled until file is loaded
			description='Please load an input file first'
		)

		self.sample_mapping_mode_container = pn.Column(
			self.sample_mapping_select,
			self.generate_samplemap_button,
			align='start'  # Align to the left
		)

		self.samplemap_fileupload = pn.widgets.FileInput(
			accept='.tsv,.csv,.txt',
			margin=(5, 5, 10, 20),
			visible=False
		)

		self.samplemap_table = pn.widgets.Tabulator(
			layout='fit_data_fill',
			height=250,
			show_index=False,
			width=500,
			margin=(5, 5, 5, 5),
			visible=False
		)

		self.assign_cond_pairs = pn.widgets.CrossSelector(
			width=600,
			height=250,
			margin=(5, 5, 10, 20),
			name='Select condition pairs'
		)

		self.modification_type = pn.widgets.TextInput(
			name='Modification type:',
			placeholder='e.g., [Phospho (STY)] for Spectronaut',
			width=300,
			description=gui_textfields.Descriptions.tooltips['ptm_settings']
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
			description=gui_textfields.Descriptions.tooltips['filtering_options']
		)

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
			"**Select** the condition pairs you want to analyze **on the left** and **move them to the right:**",
			visible=True
		)

		self.analysis_type = pn.widgets.Select(
			name='Select Condition Analysis Type',
			options=[PAIRWISE_COMPARISON, MEDIAN_CONDITION_ANALYSIS],
			value='Select an analysis',
			description='Choose between comparing pairs of conditions or comparing each condition against a median reference'
		)


		self.medianref_message = pn.pane.Markdown(
			"Every condition will be compared against the median reference",
			visible=False,
		)

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

		self.condition_progress_panel = pn.Column(
			pn.pane.Markdown("### Analysis Progress", margin=(0,0,10,0)),
			pn.pane.Markdown(
				"As soon as a condition pair is finished, you can inspect the results in the plots tabs.",
				margin=(0,0,10,0)
			),
			sizing_mode='stretch_width'
		)

		# Watchers
		self.sample_mapping_select.param.watch(self._toggle_sample_mapping_mode, 'value')
		self.path_analysis_file.param.watch(
			self._activate_after_analysis_file_upload, 'value'
		)
		self.samplemap_fileupload.param.watch(self._update_samplemap_table, 'value')
		self.samplemap_table.param.watch(self._add_conditions_for_assignment, 'value')
		self.minrep_either.param.watch(self._update_minrep_both, 'value')
		self.run_pipeline_button.param.watch(self._run_pipeline, 'clicks')
		self.analysis_type.param.watch(self._update_analysis_type_visibility, 'value')
		self.filtering_options.param.watch(self._toggle_filtering_options, 'value')
		self.path_output_folder.param.watch(self._update_results_dir, 'value')
		self.path_analysis_file.param.watch(self._update_analysis_file, 'value')
		self.samplemap_fileupload.param.watch(self._update_samplemap, 'value')
		self.generate_samplemap_button.on_click(self._generate_samplemap)
		# Add watcher for condition pairs selection
		self.assign_cond_pairs.param.watch(self._update_run_button_state, 'value')
		self.analysis_type.param.watch(self._update_run_button_state, 'value')


	def create(self):
		"""
		Build and return the main layout for the pipeline widget.
		"""

		ptm_section = pn.Row(
			pn.Column(self.modification_type, self.organism)
		)

		filtering_section = pn.Row(
			pn.Column(self.filtering_options, self.minrep_either)
		)

		advanced_settings_card = pn.Card(
			pn.Column(
				"### Threshold Settings",
				self.min_num_ions,
				self.minpep,
				self.cluster_threshold_pval,
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
			self.template_success_message,
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
			ptm_section,
			title='PTM Settings',
			collapsed=True,
			margin=(5, 5, 5, 5),
			sizing_mode='fixed',
			width=400
		)

		main_col = pn.Column(
			"### Input Files",
			self.path_analysis_file,
			self.path_output_folder,
			"### Samples and Conditions",
			samples_conditions_layout,
			self.analysis_type,
			condition_comparison_layout,
			"### Basic Settings",
			filtering_section,
			ptm_settings_card,
			advanced_settings_card,
			"### Pipeline Controls",
			pn.Row(
				self.run_pipeline_button,
				self.run_pipeline_progress,
				sizing_mode='stretch_width'
			),
			self.run_pipeline_error,
			self.condition_progress_panel,
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
		"""
		Run the alphaquant pipeline when the button is clicked.
		"""
		if self.analysis_type.value == 'Select an analysis':
			self.run_pipeline_error.object = "Please select an analysis type before running the pipeline."
			self.run_pipeline_error.visible = True
			return

		# Only check for condition pairs if doing pairwise comparison
		if self.analysis_type.value == 'Pairwise Comparison' and not self.assign_cond_pairs.value:
			self.run_pipeline_error.object = "Please select the condition pairs you want to analyze in the cross-selector above."
			self.run_pipeline_error.visible = True
			return

		# Print debug info about samplemap
		print("\n=== Pipeline Run Samplemap Check ===")
		print(f"Samplemap status: {'Present with ' + str(len(self.samplemap_table.value)) + ' rows' if self.samplemap_table.value is not None else 'None'}")
		if self.samplemap_table.value is None:
			self.run_pipeline_error.object = "Missing sample mapping. Please generate or upload a sample mapping before running the pipeline."
			self.run_pipeline_error.visible = True
			return

		print("\n=== Starting Pipeline Run ===")
		self.run_pipeline_progress.active = True
		self.run_pipeline_error.visible = False
		self.console_output.value = "Starting pipeline...\n"

		# Save samplemap to results directory
		try:
			if self.samplemap_table.value is not None and self.path_output_folder.value:
				# Create results directory if it doesn't exist
				os.makedirs(self.path_output_folder.value, exist_ok=True)

				samplemap_path = os.path.join(self.path_output_folder.value, 'samplemap.tsv')
				self.samplemap_table.value.to_csv(samplemap_path, sep='\t', index=False)
				self.console_output.value += f"Saved sample mapping to: {samplemap_path}\n"
		except Exception as e:
			self.console_output.value += f"Warning: Could not save sample mapping file: {str(e)}\n"

		# Update progress panel with selected condition pairs
		self._update_progress_panel(self.assign_cond_pairs.value or [])

		# Show overall status
		self._overall_status.visible = True

		# Set initial status for all pairs to pending
		for pair, progress_dict in self._condition_progress.items():
			print(f"Setting initial status for {pair}")
			progress_dict['pending'].visible = True
			progress_dict['running'].visible = False
			progress_dict['complete'].visible = False
			progress_dict['spinner'].visible = False

		print("Starting progress monitor...")
		try:
			# Stop any existing monitor
			if self._progress_monitor is not None:
				print("Stopping existing progress monitor...")
				self._progress_monitor.stop()

			# Create new monitor with explicit start
			self._progress_monitor = pn.state.add_periodic_callback(
				callback=self._check_condition_progress,
				period=1000,  # Check every second
				start=True
			)
			print(f"Progress monitor started successfully: {self._progress_monitor}")

			# Force an immediate check
			self._check_condition_progress()

		except Exception as e:
			print(f"Error setting up progress monitor: {str(e)}")
			import traceback
			traceback.print_exc()

		try:
			# Get condition combinations based on analysis type
			is_median_analysis = self.analysis_type.value == MEDIAN_CONDITION_ANALYSIS

			# Ensure the samplemap is saved to a file in the results directory
			samplemap_path = None
			if self.samplemap_table.value is not None and self.path_output_folder.value is not None:
				# Create the results directory if it doesn't exist
				os.makedirs(self.path_output_folder.value, exist_ok=True)

				# Save the samplemap to a file in the results directory
				samplemap_path = os.path.join(self.path_output_folder.value, "samplemap.tsv")
				self.samplemap_table.value.to_csv(samplemap_path, sep='\t', index=False)
				print(f"Saved samplemap to file: {samplemap_path}")

			if is_median_analysis:
				# For median analysis, set condition pairs to None
				cond_combinations = None
				print(f"Running median analysis with samplemap containing {len(self.samplemap_table.value)} rows")
				print(f"Samplemap columns: {', '.join(self.samplemap_table.value.columns)}")
			else:
				# For pairwise comparison, get the selected pairs
				cond_combinations = [
					tuple(pair.split(aq_variables.CONDITION_PAIR_SEPARATOR))
					for pair in self.assign_cond_pairs.value
				]

			# Log samplemap status right before passing to pipeline
			print(f"Samplemap right before pipeline run: {'Present with ' + str(len(self.samplemap_table.value)) + ' rows' if self.samplemap_table.value is not None else 'None'}")

			# Collect all configuration parameters
			pipeline_params = {
				'input_file': self.path_analysis_file.value,
				'samplemap_df': self.samplemap_table.value,
				'samplemap_file': samplemap_path,  # Pass the samplemap path to samplemap_file parameter
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
				'multicond_median_analysis': is_median_analysis,
			}

			# Log key parameters
			print("\n=== Pipeline Parameters ===")
			for key, value in pipeline_params.items():
				if key == 'samplemap_df':
					print(f"samplemap_df: {'DataFrame with ' + str(len(value)) + ' rows' if value is not None else 'None'}")
				elif key == 'samplemap_file':
					print(f"samplemap_file: {value}")
				elif key == 'condpairs_list':
					print(f"condpairs_list: {value}")
				else:
					print(f"{key}: {value}")

			# Run the pipeline
			diffmgr.run_pipeline(**pipeline_params)

			self._logger.info("Pipeline completed successfully!")

		except Exception as e:
			error_message = f"Error running pipeline: {e}"
			self.run_pipeline_error.object = error_message
			self.run_pipeline_error.visible = True
			self._logger.error(error_message)
			self._overall_status.object = "❌ Analysis Failed"
			self._overall_status.styles = {'color': 'red', 'font-weight': 'bold'}

		finally:
			print("\n=== Pipeline Run Complete ===")
			# Update overall status on completion
			if not self.run_pipeline_error.visible:
				self._overall_status.object = "✅ Analysis Complete"
				self._overall_status.styles = {'color': 'green', 'font-weight': 'bold'}
				self._overall_status.param.trigger('object')

			# Let the progress monitor run for a bit longer to catch the final status
			print("Waiting for final progress checks...")
			#time.sleep(2)  # Give time for final file checks

			# Stop progress monitoring with debug info
			if self._progress_monitor is not None:
				print("Stopping progress monitor...")
				try:
					# Force one final check before stopping
					self._check_condition_progress()
					time.sleep(1)  # Give UI time to update
					self._progress_monitor.stop()
					print("Progress monitor stopped successfully")
				except Exception as e:
					print(f"Error stopping progress monitor: {str(e)}")

			self.trigger_dependency()
			self.run_pipeline_progress.active = False

		# Show/hide components based on selected analysis type
		if self.analysis_type.value == MEDIAN_CONDITION_ANALYSIS:
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
			self.samplemap_table.visible = False
			# Hide the generate button when in upload mode
			self.generate_samplemap_button.visible = False
			# Hide the template success message when switching to upload mode
			self.template_success_message.visible = False
		else:  # 'Generate new sample to condition map'
			self.samplemap_fileupload.visible = False
			self.samplemap_table.visible = False  # Only show after button click
			# Show the generate button when in generate mode
			self.generate_samplemap_button.visible = True
			# Update button state based on whether we have an analysis file
			self._update_generate_button_state()
			# Also hide the template success message when switching back
			self.template_success_message.visible = False

	def _activate_after_analysis_file_upload(self, event):
		"""Handle analysis file upload."""
		if event.new:
			self._set_default_output_folder()
			self.path_output_folder.disabled = False
			self.run_pipeline_button.disabled = False
		else:
			self.run_pipeline_button.disabled = True

		# Update generate button state based on file presence
		self._update_generate_button_state()

	def _update_generate_button_state(self):
		"""Update generate button state based on current conditions."""
		has_file = bool(self.path_analysis_file.value)
		is_generate_mode = self.sample_mapping_select.value == 'Generate new sample to condition map'

		# Only enable the button if we have a file and are in generate mode
		self.generate_samplemap_button.disabled = not (has_file and is_generate_mode)

		# Update description based on state
		if not has_file:
			self.generate_samplemap_button.description = 'Please load an input file first'
		else:
			self.generate_samplemap_button.description = 'Generate sample mapping'

	def _set_default_output_folder(self):
		"""Set default output folder based on analysis file path."""
		print("\n=== Setting Default Output Folder ===")
		if self.path_analysis_file.value:
			base_path = os.path.dirname(self.path_analysis_file.value)
			output_path = os.path.join(base_path, 'results')
			print(f"Setting output path to: {output_path}")

			print("Updating path_output_folder widget...")
			self.path_output_folder.value = output_path

			print("Updating state...")
			self.state.results_dir = output_path
			print("Notifying subscribers...")
			self.state.notify_subscribers('results_dir')
			print("=== Finished Setting Default Output Folder ===\n")

	def _import_sample_names(self):
		if self.path_analysis_file.value:
			try:
				# Show loading indicator
				self.loading_samples_indicator.visible = True
				self.loading_samples_message.visible = True

				input_file = self.path_analysis_file.value
				_, config_dict, sep = config_dict_loader.get_input_type_and_config_dict(input_file)
				if config_dict["format"] == "longtable":
					sample_column = config_dict["sample_ID"]
					sample_names = set()
					for chunk in pd.read_csv(input_file, sep=sep, usecols=[sample_column], chunksize=400000):
						sample_names.update(chunk[sample_column].unique())
					self.sample_names = sample_names
				elif config_dict["format"] == "widetable":
					# Read the headers first to identify sample columns
					headers = pd.read_csv(input_file, sep=sep, nrows=0).columns.tolist()

					quant_pre_or_suffix = config_dict.get("quant_pre_or_suffix")
					# Filter headers to find those with the prefix or suffix
					sample_columns = [
						col for col in headers if (
							col.startswith(quant_pre_or_suffix) or
							col.endswith(quant_pre_or_suffix)
						)
					]
					self.sample_names = set([col.replace(quant_pre_or_suffix, '') for col in sample_columns])
				else:
					print("ERROR: Could not idenfity sample names in the input file.")
					self.run_pipeline_error.object = "Could not idenfity sample names . Please check your input file."
					self.run_pipeline_error.visible = True

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
		"""
		Whenever the samplemap table is updated or conditions are selected,
		update the progress panel.
		"""
		if self.samplemap_table.value is None:
			return

		df = self.samplemap_table.value
		if 'condition' in df.columns:
			unique_condit = df['condition'].dropna().unique()
			comb_condit = [
				aq_variables.CONDITION_PAIR_SEPARATOR.join(comb)
				for comb in itertools.permutations(unique_condit, 2)
			]
			self.assign_cond_pairs.options = comb_condit

			# Only update progress panel if pipeline is running
			if self._progress_monitor is not None:
				self._update_progress_panel(self.assign_cond_pairs.value or [])

	def _update_progress_panel(self, condition_pairs):
		"""Create or update progress indicators for each condition pair."""
		# Only show progress for selected condition pairs
		selected_pairs = self.assign_cond_pairs.value or []

		# Clear existing progress indicators
		self._condition_progress.clear()
		progress_items = []

		# Add overall progress indicator
		self._overall_status = pn.pane.Markdown(
			"⏳ Analysis Pipeline Running...",
			styles={
				'color': '#007bff',  # Bootstrap primary blue
				'font-weight': 'bold',
				'margin-bottom': '10px'
			},
			visible=False
		)
		progress_items.append(self._overall_status)

		for pair in selected_pairs:
			# Create status indicator with both spinner and text
			status_row = pn.Row(
				pn.indicators.LoadingSpinner(value=False, color='primary', width=20, height=20),
				pn.pane.Markdown('⏳', styles={'font-size': '20px'}, visible=False),  # Pending
				pn.pane.Markdown('🔄', styles={'font-size': '20px'}, visible=False),  # Running
				pn.pane.Markdown('✅', styles={'color': 'green', 'font-size': '20px'}, visible=False),  # Complete
				pn.pane.Markdown(f"**{pair}**"),
				margin=(5, 5, 5, 10)
			)

			self._condition_progress[pair] = {
				'row': status_row,
				'spinner': status_row[0],
				'pending': status_row[1],
				'running': status_row[2],
				'complete': status_row[3]
			}
			progress_items.append(status_row)

		# Create a card with the progress indicators
		progress_card = pn.Card(
			pn.Column(
				pn.pane.Markdown("### Analysis Progress"),
				pn.pane.Markdown(
					"As soon as a condition pair is finished, you can inspect the results in the plots tabs.",
					margin=(0, 0, 10, 0)
				),
				*progress_items,
				margin=(10, 10, 10, 10)
			),
			title='Progress',
			collapsed=False,
			margin=(5, 5, 5, 5),
			styles={
				'background': '#f8f9fa',
				'border': '1px solid #dee2e6',
				'border-radius': '5px',
				'box-shadow': '0 1px 3px rgba(0,0,0,0.12)'
			}
		)

		# Update the progress panel
		self.condition_progress_panel.clear()
		self.condition_progress_panel.append(progress_card)

	def _check_condition_progress(self):
		"""Check progress of condition pairs by looking for result files."""
		try:
			print("\n=== Checking Condition Progress ===")

			if self.path_output_folder.value is None:
				print("No output folder value, returning early")
				return

			if not self._condition_progress:
				print("No condition progress dictionary, returning early")
				return

			print(f"Checking output folder: {self.path_output_folder.value}")
			print(f"Number of condition pairs to check: {len(self._condition_progress)}")

			all_complete = True

			for pair, progress_dict in self._condition_progress.items():
				print(f"\nChecking pair: {pair}")
				cond1, cond2 = pair.split('_VS_')
				result_file = os.path.join(
					self.path_output_folder.value,
					f"{cond1}_VS_{cond2}.results.tsv"
				)
				print(f"Looking for file: {result_file}")

				file_exists = os.path.exists(result_file)
				print(f"File exists: {file_exists}")

				if file_exists:
					print(f"Found result file for {pair}")
					progress_dict['pending'].visible = False
					progress_dict['running'].visible = False
					progress_dict['complete'].visible = True
					progress_dict['spinner'].visible = False
				else:
					print(f"Result file not found for {pair}")
					progress_dict['pending'].visible = False
					progress_dict['running'].visible = True
					progress_dict['complete'].visible = False
					progress_dict['spinner'].visible = True
					all_complete = False

			print(f"\nAll conditions complete? {all_complete}")
			print(f"Current overall status: {self._overall_status.object}")

			if all_complete and self._overall_status.object == "⏳ Analysis Pipeline Running...":
				print("Updating overall status to complete")
				self._overall_status.object = "✅ Analysis Complete"
				self._overall_status.styles = {'color': 'green', 'font-weight': 'bold'}
				self._overall_status.param.trigger('object')

			print("=== Progress Check Complete ===\n")

		except Exception as e:
			print(f"Error in progress check: {str(e)}")
			import traceback
			traceback.print_exc()

	def _update_minrep_both(self, *events):
		"""Set minrep_both to 0 when minrep_either is changed."""
		self.minrep_both.value = 0

	def _update_results_dir(self, event):
		"""Update central state with new results directory."""
		print("\n=== Updating Results Directory ===")
		print(f"Event type: {type(event)}")
		print(f"Event value: {event}")
		if hasattr(event, 'new'):
			print(f"Event.new value: {event.new}")
		self.state.results_dir = event.new
		print("Notifying subscribers...")
		self.state.notify_subscribers('results_dir')
		print("=== Finished Updating Results Directory ===\n")

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
				print(f"Error reading sample map: {str(e)}")
				self.run_pipeline_error.object = f"Error reading sample map: {str(e)}"
				self.run_pipeline_error.visible = True

	def _generate_samplemap(self, event):
		"""Handle the generate button click event."""
		# Show loading indicators
		self.loading_samples_indicator.visible = True
		self.loading_samples_message.visible = True
		self.template_success_message.visible = False  # Hide any previous success message

		try:
			self._import_sample_names()
			self._init_samplemap_df_template()
			self.samplemap_table.visible = True

			if self.path_output_folder.value:  # Only save if output folder is set
				os.makedirs(self.path_output_folder.value, exist_ok=True)
				template_path = os.path.join(self.path_output_folder.value, 'samplemap_template.tsv')
				self.samplemap_table.value.to_csv(
					template_path,
					sep="\t",
					index=None
				)
				print("wrote samplemap template to disk")

				# Show success message with the actual path
				self.template_success_message.object = f"""Template has been generated. Please fill out the condition column in the table below.\nThe template has also been saved to
				<code>{template_path}</code>\nif you prefer to edit it with Excel or other applications."""
				self.template_success_message.visible = True
		except Exception as e:
			print(f"Error generating sample map: {str(e)}")
			self.run_pipeline_error.object = f"Error generating sample map: {str(e)}"
			self.run_pipeline_error.visible = True
		finally:
			# Hide loading indicators when done
			self.loading_samples_indicator.visible = False
			self.loading_samples_message.visible = False

	def natural_sort(self, l):
		"""
		Sort a list in a way that numbers are sorted numerically rather than alphabetically.
		"""
		convert = lambda text: int(text) if text.isdigit() else text.lower()
		alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
		return sorted(l, key=alphanum_key)

	def _update_analysis_type_visibility(self, change=None):
		"""
		Update the visibility of UI elements based on the selected analysis type.
		"""
		# Handle both dictionary-like objects and param.parameterized.Event objects
		if change is None:
			pass
		elif hasattr(change, 'get'):  # Dictionary-like object
			if change.get('name') != 'value':
				return
		elif hasattr(change, 'name'):  # param.parameterized.Event object
			if change.name != 'value':
				return
		else:
			return

		analysis_type = self.analysis_type.value

		# Show/hide condition pairs selector and related elements based on analysis type
		if analysis_type == 'Pairwise Comparison':
			self.assign_cond_pairs.visible = True
			self.condition_comparison_header.visible = True
			self.condition_comparison_instructions.visible = True
			self.medianref_message.visible = False
		elif analysis_type == MEDIAN_CONDITION_ANALYSIS:
			self.assign_cond_pairs.visible = False
			self.condition_comparison_header.visible = False
			self.condition_comparison_instructions.visible = False
			self.medianref_message.visible = True
		else:
			self.assign_cond_pairs.visible = False
			self.condition_comparison_header.visible = False
			self.condition_comparison_instructions.visible = False
			self.medianref_message.visible = False

		# Update the run button state
		self._update_run_button_state()

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
		self._log_stream.seek(0)
		new_logs = self._log_stream.read()

		# Update the console widget
		if new_logs:
			current_text = self.console_output.value or ""
			self.console_output.value = current_text + new_logs

			# Clear the StringIO buffer
			self._log_stream.truncate(0)
			self._log_stream.seek(0)

	def __del__(self):
		"""Clean up logging handler when the widget is destroyed."""
		if self._stream_handler is not None:
			self._logger.removeHandler(self._stream_handler)
			self._stream_handler.close()

	def _update_run_button_state(self, event=None):
		"""
		Update the run button state based on condition pairs selection and analysis type.
		"""
		# For Median Condition Analysis, we don't need condition pairs
		if self.analysis_type.value == MEDIAN_CONDITION_ANALYSIS:
			self.run_pipeline_button.disabled = False
			self.run_pipeline_button.description = 'Run pipeline'
		# For Pairwise Comparison, we need at least one condition pair
		elif self.analysis_type.value == 'Pairwise Comparison':
			has_pairs = bool(self.assign_cond_pairs.value)
			self.run_pipeline_button.disabled = not has_pairs
			if not has_pairs:
				self.run_pipeline_button.description = 'Please select condition pairs'
			else:
				self.run_pipeline_button.description = 'Run pipeline'
		# If no analysis type is selected
		else:
			self.run_pipeline_button.disabled = True
			self.run_pipeline_button.description = 'Please select an analysis type'

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
	# Create state manager first
	state_manager = gui.DashboardState()

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

	# Create pipeline instance with state manager
	pipeline = RunPipeline(state=state_manager)
	pipeline_layout = pipeline.create()

	# Create plotting tabs with state manager
	plotting_tab = dashboad_parts_plots_basic.PlottingTab(state=state_manager)
	proteoform_tab = dashboad_parts_plots_proteoforms.ProteoformPlottingTab(state=state_manager)

	# Register subscribers
	state_manager.register_subscriber(plotting_tab)
	state_manager.register_subscriber(proteoform_tab)

	# Create tabs
	all_tabs = pn.Tabs(
		('Pipeline', pipeline_layout),
		('Single Comparison', plotting_tab.panel()),
		('Plotting', proteoform_tab.panel()),
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
