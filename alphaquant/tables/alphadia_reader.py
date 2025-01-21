import numpy as np
import pandas as pd
from typing import Dict
import alphaquant.utils.utils as aq_utils
import os

import alphaquant.config.config as aqconfig
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

class AlphaDIAFragTableProcessor:
	"""
	A class for processing DIA fragment tables and generating quantification data.

	This class handles the processing of fragment matrices from DIA mass spectrometry data,
	including precursor mapping, quantification ID generation, and ML info table creation.
	"""

	def __init__(self, fragment_matrix_file: str):
		"""
		Initialize the AlphaDIAFragTableProcessor.

		Args:
			fragment_matrix_file (str): Path to the fragment matrix file
		"""

		self.ml_info_file = aq_utils.get_progress_folder_filename(fragment_matrix_file, ".ml_info_table.tsv")
		self.input_file_reformat = aq_utils.get_progress_folder_filename(fragment_matrix_file, ".alphadia_fragion.aq_reformat.tsv", remove_extension=False)

		precursor_file = os.path.join(os.path.dirname(fragment_matrix_file), "precursors.tsv")
		self._precursor_df = pd.read_csv(precursor_file, sep="\t")
		self._precursor2quantID = self._precursor2quantid()

		if not os.path.exists(self.ml_info_file):
			LOGGER.info(f"Creating ML info file")
			self.ml_info_df = self._define_ml_info_table()
			self.ml_info_df.to_csv(self.ml_info_file, sep="\t", index=False)
		else:
			LOGGER.info(f"ML info file already exists at {self.ml_info_file}")

		if not os.path.exists(self.input_file_reformat):
			LOGGER.info(f"Creating processed fragment matrix")
			self._fragment_matrix_df = pd.read_csv(fragment_matrix_file, sep="\t")
			self.processed_fragment_matrix = self._process_fragment_matrix()
			self.processed_fragment_matrix.to_csv(self.input_file_reformat, sep="\t", index=False)
		else:
			LOGGER.info(f"Processed fragment matrix already exists at {self.input_file_reformat}")

	def _precursor2quantid(self) -> Dict[str, str]:
		"""
		Generate quantification IDs by combining sequence, modification, and charge information.
		Updates the precursor_df with new quant_id column and creates a mapping dictionary.

		Returns:
			Dict[str, str]: Mapping of mod_seq_charge_hash to quant_id
		"""
		# Clean modification data
		mods_clean = self._precursor_df['mods'].fillna('')
		mod_sites_clean = self._precursor_df['mod_sites'].fillna('')

		# Combine modifications and sites
		mod_combined = (mods_clean + '_' + mod_sites_clean).where(
			(mods_clean != '') & (mod_sites_clean != ''),
			''
		)

		# Generate quant_id
		self._precursor_df["quant_id"] = (
			'SEQ_' + self._precursor_df['sequence'] +
			'_MOD_' + mod_combined +
			'_CHARGE_' + self._precursor_df['charge'].astype(str)
		)

		# Create mapping dictionary using mod_seq_charge_hash
		self._precursor2quantID = dict(zip(
			self._precursor_df['mod_seq_charge_hash'],
			self._precursor_df['quant_id']
		))

		return self._precursor2quantID

	def _define_ml_info_table(self) -> pd.DataFrame:
		"""
		Create a machine learning information table from the precursor data.

		Returns:
			pd.DataFrame: Processed DataFrame containing only numeric columns and essential identifiers
		"""
		if "quant_id" not in self._precursor_df.columns:
			raise ValueError("Must call _precursor2quantid() before creating ML info table")

		ml_info_df = (
			self._precursor_df
			.set_index(["quant_id", "run"])
			.select_dtypes(include=[np.number])
			.reset_index()
			.rename(columns={"run": "sample_ID"})
		)
		return ml_info_df

	def _process_fragment_matrix(self) -> pd.DataFrame:
		"""
		Process the fragment matrix by mapping quantification IDs and reformatting columns.

		Returns:
			pd.DataFrame: Processed fragment matrix in the desired format
		"""
		if self._precursor2quantID is None:
			raise ValueError("Must call _precursor2quantid() before processing fragment matrix")

		# Map quant_ids to fragment matrix using mod_seq_charge_hash
		processed_fragment_matrix_df = self._fragment_matrix_df.copy()
		processed_fragment_matrix_df["quant_id"] = (
			processed_fragment_matrix_df['mod_seq_charge_hash'].map(self._precursor2quantID) +
			"_FRGION_" + processed_fragment_matrix_df["ion"].astype(str)
		)

		# Format for AQ
		columns_to_drop = ["ion", "mod_seq_hash", "mod_seq_charge_hash", "precursor_idx"]
		processed_fragment_matrix_df = processed_fragment_matrix_df.drop(columns=columns_to_drop).rename(columns={"pg": "protein"})

		return processed_fragment_matrix_df
