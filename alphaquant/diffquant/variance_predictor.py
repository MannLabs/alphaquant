"""Variance-predictor scoring for background-distribution ion ordering.

Loads quality-metric columns from the ml_info_table (at precursor level)
and fits a linear model predicting observed per-ion variance from these
metrics.  The predicted variance is used to sort ions so that those with
similar expected variance are grouped together, producing tighter
empirical background distributions.

The linear-regression approach automatically handles columns that
correlate with variance in different directions (positive or negative
coefficients) and weights each column by its predictive power.
"""

import re
import logging

import numpy as np
import pandas as pd

import alphaquant.utils.reader_utils as aq_reader_utils

LOGGER = logging.getLogger(__name__)

_ION_SPLIT_PAT = re.compile(r"(_FRGION_|_MS1ISOTOPES_|_PRECURSOR_)")


def load_variance_predictor_scores(ml_info_file, samples_used,
                                   variance_predictor_cols, ion_index,
                                   ion_variance, ion_median_intensity=None):
    """Fit a linear model predicting observed ion variance from quality
    metrics and return the predicted scores for sorting.

    The ml_info_table is at precursor level (SEQ_MOD_CHARGE).  Fragment-level
    ion ids are mapped to their parent precursor via ``_ION_SPLIT_PAT``.

    A simple OLS regression is fitted::

        ion_variance ~ median_intensity + col1 + col2 + ...

    Median intensity is always included as a built-in predictor (when
    provided) because it is the strongest single correlate of ion-level
    variance.  The quality-metric columns from the config refine this
    baseline.  The predicted values serve as the combined score.

    Args:
        ml_info_file (str): Path to the ml_info_table TSV.
        samples_used (list[str] | None): Sample IDs to keep; None keeps all.
        variance_predictor_cols (list[str]): Column names in the ml_info_table
            to use as predictors.
        ion_index (pd.Index): Ion identifiers to score.
        ion_variance (pd.Series): Observed per-ion variance (indexed by ion id),
            used as the regression target.
        ion_median_intensity (pd.Series | None): Per-ion pooled median
            intensity.  When provided, prepended as a built-in predictor.

    Returns:
        dict[str, float] | None: Mapping from ion id to predicted variance
            score, or None if the required columns were not found or the
            regression could not be fitted.
    """
    if not variance_predictor_cols:
        return None

    try:
        usecols = ["quant_id", "sample_ID"] + list(variance_predictor_cols)
        ml_df = aq_reader_utils.read_file(ml_info_file, sep="\t", usecols=usecols)
    except (ValueError, KeyError):
        LOGGER.warning(
            "Could not load variance predictor columns %s from %s. "
            "Falling back to intensity-only background sorting.",
            variance_predictor_cols, ml_info_file,
        )
        return None

    if samples_used is not None:
        ml_df = ml_df[ml_df["sample_ID"].isin(samples_used)]
    ml_df = ml_df.drop(columns=["sample_ID"])

    for col in variance_predictor_cols:
        ml_df[col] = pd.to_numeric(ml_df[col], errors="coerce")

    available_cols = [c for c in variance_predictor_cols if c in ml_df.columns]
    if not available_cols:
        LOGGER.warning("None of the variance predictor columns found. Falling back.")
        return None

    prec_features = ml_df.groupby("quant_id")[available_cols].median()

    ion_strings = np.array([str(x) for x in ion_index])
    precursor_ids = np.array([_ION_SPLIT_PAT.split(s)[0] for s in ion_strings])

    scores = _fit_and_predict(prec_features, available_cols,
                              precursor_ids, ion_index, ion_variance,
                              ion_median_intensity)
    if scores is None:
        return None

    all_cols = (["median_intensity"] if ion_median_intensity is not None
                else []) + available_cols
    ion2varscore = dict(zip(ion_index, scores))
    LOGGER.info(
        "Variance predictor scores computed for %d/%d ions using columns %s",
        np.isfinite(scores).sum(), len(ion_index), all_cols,
    )
    return ion2varscore


def _fit_and_predict(prec_features, available_cols, precursor_ids,
                     ion_index, ion_variance,
                     ion_median_intensity=None):
    """Build feature matrix, fit OLS on observed variance, return predictions.

    Median intensity (when provided) is always prepended as the first
    predictor column because it is typically the strongest correlate of
    ion-level variance.  The quality-metric columns refine the prediction.

    Ions with missing features or missing variance are excluded from fitting
    but receive the median predicted score as fallback.

    Args:
        prec_features (pd.DataFrame): Precursor-level feature medians
            (index = quant_id, columns = available_cols).
        available_cols (list[str]): Column names to use as predictors.
        precursor_ids (np.ndarray): Precursor id for each ion in ion_index.
        ion_index (pd.Index): Ion identifiers.
        ion_variance (pd.Series): Observed per-ion variance.
        ion_median_intensity (pd.Series | None): Per-ion pooled median
            intensity.  When provided, used as a built-in first predictor.

    Returns:
        np.ndarray of predicted scores (length = len(ion_index)), or None
        if the regression cannot be fitted.
    """
    n_ions = len(ion_index)
    n_cols = len(available_cols)

    X = np.full((n_ions, n_cols), np.nan)
    for j, col in enumerate(available_cols):
        col_vals = prec_features[col]
        X[:, j] = [col_vals.get(pid, np.nan) for pid in precursor_ids]

    if ion_median_intensity is not None:
        intensity_col = np.array(
            [ion_median_intensity.get(ion, np.nan) for ion in ion_index],
            dtype=float,
        ).reshape(-1, 1)
        X = np.column_stack([intensity_col, X])
        all_col_names = ["median_intensity"] + list(available_cols)
    else:
        all_col_names = list(available_cols)

    n_features = X.shape[1]

    y = np.array([ion_variance.get(ion, np.nan) for ion in ion_index],
                 dtype=float)

    valid = np.isfinite(X).all(axis=1) & np.isfinite(y)
    if valid.sum() < max(10, n_features + 1):
        LOGGER.warning(
            "Too few valid ions (%d) for variance predictor regression. "
            "Falling back to intensity-only sorting.", valid.sum()
        )
        return None

    X_fit = X[valid]
    y_fit = y[valid]

    # Standardise features for numerical stability
    X_mean = X_fit.mean(axis=0)
    X_std = X_fit.std(axis=0)
    X_std[X_std < 1e-15] = 1.0
    X_fit_z = (X_fit - X_mean) / X_std

    # OLS with intercept: y = X_z @ beta + intercept
    X_design = np.column_stack([np.ones(X_fit_z.shape[0]), X_fit_z])
    try:
        beta, _, _, _ = np.linalg.lstsq(X_design, y_fit, rcond=None)
    except np.linalg.LinAlgError:
        LOGGER.warning("OLS fit failed for variance predictor. Falling back.")
        return None

    LOGGER.info(
        "Variance predictor coefficients (standardised): %s",
        dict(zip(all_col_names, beta[1:])),
    )

    # Predict for all ions (including those excluded from fit)
    X_all_z = (X - X_mean) / X_std
    X_all_design = np.column_stack([np.ones(n_ions), X_all_z])
    predicted = X_all_design @ beta

    # Fallback for ions with missing features
    finite_mask = np.isfinite(predicted)
    if finite_mask.any():
        median_pred = np.median(predicted[finite_mask])
    else:
        median_pred = 0.0
    predicted = np.where(finite_mask, predicted, median_pred)

    return predicted
