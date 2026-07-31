"""Peptide-level intensity aggregation via directLFQ-style trace alignment.

Groups base ions (fragment ions, MS1 isotopes, precursors) by their parent
peptide (mod_seq_charge level), aligns the ion intensity traces across
samples, and derives a single peptide intensity per sample.

Algorithm per peptide group:
  1. Collect all ion traces (log2 intensities x samples).
  2. Align traces using hierarchical agglomerative matching (same
     algorithm as within-condition sample normalization).
  3. Take the per-sample median of aligned traces.
  4. Scale so that total linear intensity is preserved.
"""

import re
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd

from alphaquant.norm.normalization import get_normfacts_withincond, apply_sampleshifts

import alphaquant.config.config as aqconfig
import logging

aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

_PEPTIDE_RE = re.compile(r"(SEQ.*MOD.*CHARGE.*)(_FRG.*|_MS1.*|_PREC.*)")


def _get_peptide_id(ion_name):
    """Extract the peptide (mod_seq_charge) identifier from a base ion name."""
    m = _PEPTIDE_RE.match(ion_name)
    if m:
        return m.group(1)
    return ion_name


def _align_traces(ion_matrix):
    """Align ion traces using hierarchical agglomerative matching.

    Each row is one ion's intensity profile across samples (log2).
    Rows are shifted vertically so that they overlap optimally.

    Args:
        ion_matrix: (n_ions, n_samples) numpy array in log2 space.

    Returns:
        Aligned copy of the matrix (same shape).
    """
    aligned = ion_matrix.copy()
    if aligned.shape[0] <= 1:
        return aligned
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        idx2shift = get_normfacts_withincond(aligned)
    apply_sampleshifts(aligned, idx2shift)
    return aligned


def _derive_peptide_intensity(ion_matrix_original, ion_matrix_aligned, min_nonan=1):
    """Derive a single peptide intensity profile from aligned ion traces.

    Takes the per-sample median of the aligned traces and scales the result
    so that the total linear intensity is preserved.

    Args:
        ion_matrix_original: (n_ions, n_samples) log2 values before alignment.
        ion_matrix_aligned:  (n_ions, n_samples) log2 values after alignment.
        min_nonan: minimum non-NaN ions per sample to produce a value.

    Returns:
        1-D array of length n_samples (log2 peptide intensities), or None
        if the result would be all-NaN.
    """
    n_samples = ion_matrix_aligned.shape[1]

    total_original = np.nansum(2.0 ** ion_matrix_original)

    peptide_profile = np.full(n_samples, np.nan)
    for s in range(n_samples):
        col = ion_matrix_aligned[:, s]
        if np.sum(~np.isnan(col)) >= min_nonan:
            peptide_profile[s] = np.nanmedian(col)

    total_profile = np.nansum(2.0 ** peptide_profile)
    if total_profile == 0:
        return None

    scale = total_original / total_profile
    peptide_profile += np.log2(scale)
    return peptide_profile


def aggregate_ions_to_peptides(df_combined, pep2prot):
    """Aggregate base-ion intensities to peptide level via trace alignment.

    Both conditions should be combined into a single DataFrame before calling
    this function so that ion traces are aligned consistently across all
    samples.

    Args:
        df_combined: DataFrame with log2 intensities.
            Index = base-ion names, columns = sample names.
        pep2prot: dict mapping base-ion name -> protein name.

    Returns:
        (df_peptides, pep2prot_new):
            df_peptides: DataFrame with peptide-level log2 intensities
                (index = peptide ids, columns = sample names).
            pep2prot_new: dict mapping peptide id -> protein name.
    """
    peptide2ions = defaultdict(list)
    for ion_name in df_combined.index:
        pep_id = _get_peptide_id(ion_name)
        peptide2ions[pep_id].append(ion_name)

    LOGGER.info(
        "Aggregating %d base ions into %d peptides via trace alignment",
        len(df_combined.index),
        len(peptide2ions),
    )

    peptide_names = []
    peptide_rows = []
    pep2prot_new = {}

    for pep_id, ion_names in peptide2ions.items():
        present = [ion for ion in ion_names if ion in df_combined.index]
        if not present:
            continue

        ion_matrix = df_combined.loc[present].to_numpy(dtype=float)
        if np.all(np.isnan(ion_matrix)):
            continue

        aligned = _align_traces(ion_matrix)
        profile = _derive_peptide_intensity(ion_matrix, aligned)
        if profile is None:
            continue

        peptide_names.append(pep_id)
        peptide_rows.append(profile)

        for ion in present:
            if ion in pep2prot:
                pep2prot_new[pep_id] = pep2prot[ion]
                break

    df_peptides = pd.DataFrame(
        peptide_rows, index=peptide_names, columns=df_combined.columns,
    )

    LOGGER.info("Peptide aggregation complete: %d peptides", len(df_peptides))
    return df_peptides, pep2prot_new


def aggregate_conditions_to_peptides(df_c1, df_c2, pep2prot):
    """Convenience wrapper: aggregate two per-condition DataFrames at once.

    Combines both conditions so that traces are aligned across all samples,
    runs aggregation, then splits back into per-condition DataFrames.

    Args:
        df_c1: condition-1 DataFrame (ions x c1_samples), log2 intensities.
        df_c2: condition-2 DataFrame (ions x c2_samples), log2 intensities.
        pep2prot: dict mapping base-ion name -> protein name.

    Returns:
        (df_c1_pep, df_c2_pep, pep2prot_new)
    """
    df_combined = df_c1.join(df_c2, how="outer")
    df_peptides, pep2prot_new = aggregate_ions_to_peptides(df_combined, pep2prot)

    df_c1_pep = df_peptides[df_c1.columns]
    df_c2_pep = df_peptides[df_c2.columns]
    return df_c1_pep, df_c2_pep, pep2prot_new
