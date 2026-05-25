import numpy as np
import pandas as pd
import pytest

import alphaquant.diffquant.variance_predictor as aq_vp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

N = 20  # enough ions to pass the min-10 threshold


def _precursor_ids(n=N):
    return [f"PREC{i}_MOD{i}_{i}" for i in range(n)]


def _ion_index(precursor_ids):
    """One fragment ion per precursor — keeps precursor mapping trivial."""
    return pd.Index([f"{pid}_FRGION_y1" for pid in precursor_ids])


def _write_ml_info(tmp_path, precursor_ids, columns_dict, samples=("S1", "S2")):
    """Write a synthetic ml_info_table TSV and return its path."""
    rows = []
    for i, pid in enumerate(precursor_ids):
        for s in samples:
            row = {"quant_id": pid, "sample_ID": s}
            for col, vals in columns_dict.items():
                row[col] = vals[i]
            rows.append(row)
    df = pd.DataFrame(rows)
    path = tmp_path / "ml_info.tsv"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


def _make_ion_variance(ion_index, values):
    """Build an ion_variance Series aligned with ion_index."""
    return pd.Series(values, index=ion_index)


def _make_ion_median_intensity(ion_index, values):
    """Build an ion_median_intensity Series aligned with ion_index."""
    return pd.Series(values, index=ion_index)


def _spearman(a, b):
    from scipy.stats import spearmanr
    r, _ = spearmanr(a, b)
    return r


# ---------------------------------------------------------------------------
# Basic interface tests
# ---------------------------------------------------------------------------

class TestLoadVariancePredictorScoresInterface:
    def test_returns_none_when_no_columns(self, tmp_path):
        pids = _precursor_ids()
        ion_idx = _ion_index(pids)
        ml_file = _write_ml_info(tmp_path, pids, {"col": list(range(N))})
        ion_var = _make_ion_variance(ion_idx, np.ones(len(ion_idx)))
        result = aq_vp.load_variance_predictor_scores(
            ml_file, ["S1", "S2"], variance_predictor_cols=[],
            ion_index=ion_idx, ion_variance=ion_var,
        )
        assert result is None

    def test_returns_none_when_missing_columns(self, tmp_path):
        pids = _precursor_ids()
        ion_idx = _ion_index(pids)
        ml_file = _write_ml_info(tmp_path, pids, {"col": list(range(N))})
        ion_var = _make_ion_variance(ion_idx, np.ones(len(ion_idx)))
        result = aq_vp.load_variance_predictor_scores(
            ml_file, ["S1", "S2"], variance_predictor_cols=["NoSuchCol"],
            ion_index=ion_idx, ion_variance=ion_var,
        )
        assert result is None

    def test_returns_dict_with_correct_keys(self, tmp_path):
        pids = _precursor_ids()
        ion_idx = _ion_index(pids)
        col_vals = [float(i) for i in range(N)]
        ion_var = _make_ion_variance(ion_idx, col_vals)
        ml_file = _write_ml_info(tmp_path, pids, {"score": col_vals})
        result = aq_vp.load_variance_predictor_scores(
            ml_file, ["S1", "S2"], variance_predictor_cols=["score"],
            ion_index=ion_idx, ion_variance=ion_var,
        )
        assert result is not None
        assert set(result.keys()) == set(ion_idx)

    def test_ions_from_same_precursor_get_same_score(self, tmp_path):
        pids = _precursor_ids()
        ions = []
        for pid in pids:
            ions.append(f"{pid}_FRGION_y3")
            ions.append(f"{pid}_FRGION_y4")
        ion_idx = pd.Index(ions)
        col_vals = [float(i) for i in range(N)]
        var_vals = [float(i) * 0.1 for i in range(N)]
        # Each pair of ions from the same precursor gets the same variance
        ion_var_vals = []
        for v in var_vals:
            ion_var_vals.extend([v, v])
        ion_var = _make_ion_variance(ion_idx, ion_var_vals)
        ml_file = _write_ml_info(tmp_path, pids, {"score": col_vals})
        result = aq_vp.load_variance_predictor_scores(
            ml_file, ["S1", "S2"], variance_predictor_cols=["score"],
            ion_index=ion_idx, ion_variance=ion_var,
        )
        assert result is not None
        for pid in pids:
            assert result[f"{pid}_FRGION_y3"] == result[f"{pid}_FRGION_y4"]


# ---------------------------------------------------------------------------
# Regression logic tests with controlled distributions
# ---------------------------------------------------------------------------

class TestRegressionLogic:
    """Verify the linear-regression approach correctly recovers the
    relationship between quality metrics and ion variance."""

    def _get_scores(self, tmp_path, columns_dict, variance_predictor_cols,
                    ion_var_values, ion_med_int_values=None):
        pids = _precursor_ids()
        ion_idx = _ion_index(pids)
        ml_file = _write_ml_info(tmp_path, pids, columns_dict)
        ion_var = _make_ion_variance(ion_idx, ion_var_values)
        ion_med_int = (_make_ion_median_intensity(ion_idx, ion_med_int_values)
                       if ion_med_int_values is not None else None)
        result = aq_vp.load_variance_predictor_scores(
            ml_file, ["S1", "S2"],
            variance_predictor_cols=variance_predictor_cols,
            ion_index=ion_idx, ion_variance=ion_var,
            ion_median_intensity=ion_med_int,
        )
        assert result is not None
        return np.array([result[ion] for ion in ion_idx])

    # -- single column, positively correlated with variance ----------------

    def test_single_column_positive_correlation(self, tmp_path):
        """Column that increases with variance → scores should track variance."""
        col_vals = [float(i) for i in range(N)]
        true_var = [float(i) * 0.5 for i in range(N)]  # same direction
        scores = self._get_scores(tmp_path, {"col": col_vals}, ["col"], true_var)
        rho = _spearman(true_var, scores)
        assert rho > 0.95

    # -- single column, negatively correlated with variance ----------------

    def test_single_column_negative_correlation(self, tmp_path):
        """Column that decreases with variance → negative coefficient,
        but predicted scores should still track variance."""
        col_vals = [float(N - 1 - i) for i in range(N)]  # high col = low var
        true_var = [float(i) * 0.5 for i in range(N)]
        scores = self._get_scores(tmp_path, {"col": col_vals}, ["col"], true_var)
        rho = _spearman(true_var, scores)
        assert rho > 0.95

    # -- two columns in OPPOSITE directions: the key improvement -----------

    def test_two_columns_opposite_direction_both_recovered(self, tmp_path):
        """One column goes up with variance, the other goes down.
        The regression should handle both correctly and produce scores
        that correlate strongly with the true variance."""
        true_var = [float(i) for i in range(N)]
        col_up = [float(i) * 2.0 for i in range(N)]             # positive corr
        col_down = [float(N - 1 - i) * 3.0 for i in range(N)]   # negative corr
        scores = self._get_scores(
            tmp_path, {"col_up": col_up, "col_down": col_down},
            ["col_up", "col_down"], true_var,
        )
        rho = _spearman(true_var, scores)
        assert rho > 0.95, (
            f"Expected strong correlation even with opposing columns, got rho={rho:.3f}"
        )

    # -- magnitude irrelevant (regression handles it) ----------------------

    def test_different_magnitudes(self, tmp_path):
        """Columns with very different scales should be handled by the
        standardisation in the regression."""
        true_var = [float(i) for i in range(N)]
        col_tiny = [0.001 * i for i in range(N)]
        col_huge = [1e6 * i for i in range(N)]
        scores = self._get_scores(
            tmp_path, {"tiny": col_tiny, "huge": col_huge},
            ["tiny", "huge"], true_var,
        )
        rho = _spearman(true_var, scores)
        assert rho > 0.95

    # -- noisy predictor ---------------------------------------------------

    def test_noisy_predictor_still_useful(self, tmp_path):
        """A clean + noisy column, both in the same direction: combined
        should still correlate well with true variance."""
        rng = np.random.RandomState(42)
        true_var = [float(i) for i in range(N)]
        col_clean = [float(i) for i in range(N)]
        col_noisy = [float(i) + rng.normal(0, 5) for i in range(N)]
        scores = self._get_scores(
            tmp_path, {"clean": col_clean, "noisy": col_noisy},
            ["clean", "noisy"], true_var,
        )
        rho = _spearman(true_var, scores)
        assert rho > 0.8

    # -- uninformative column gets ~zero weight ----------------------------

    def test_uninformative_column_ignored(self, tmp_path):
        """A random column uncorrelated with variance should not hurt
        when combined with a good predictor."""
        rng = np.random.RandomState(99)
        true_var = [float(i) for i in range(N)]
        col_good = [float(i) for i in range(N)]
        col_random = list(rng.randn(N))
        scores_combined = self._get_scores(
            tmp_path, {"good": col_good, "random": col_random},
            ["good", "random"], true_var,
        )
        scores_good_only = self._get_scores(
            tmp_path, {"good": col_good}, ["good"], true_var,
        )
        rho_combined = _spearman(true_var, scores_combined)
        rho_good = _spearman(true_var, scores_good_only)
        # Combined should be almost as good as good-only
        assert rho_combined > 0.9
        assert rho_combined >= rho_good - 0.1

    # -- three columns, two agree, one opposes -----------------------------

    def test_three_columns_mixed_directions(self, tmp_path):
        """Two columns positively, one negatively associated with variance.
        Regression should handle all three correctly."""
        true_var = [float(i) for i in range(N)]
        col_up1 = [float(i) for i in range(N)]
        col_up2 = [float(i) * 3.0 for i in range(N)]
        col_down = [float(N - 1 - i) for i in range(N)]
        scores = self._get_scores(
            tmp_path,
            {"up1": col_up1, "up2": col_up2, "down": col_down},
            ["up1", "up2", "down"], true_var,
        )
        rho = _spearman(true_var, scores)
        assert rho > 0.95

    # -- nonlinear relationship: regression captures the linear component --

    def test_nonlinear_variance_partially_captured(self, tmp_path):
        """Even if the true relationship is quadratic, the linear model
        should capture the monotonic trend reasonably well."""
        true_var = [float(i)**2 for i in range(N)]
        col_vals = [float(i) for i in range(N)]  # linear predictor
        scores = self._get_scores(tmp_path, {"col": col_vals}, ["col"], true_var)
        rho = _spearman(true_var, scores)
        # Spearman only cares about rank order, so linear model on monotonic
        # data should give perfect rank correlation
        assert rho > 0.95


# ---------------------------------------------------------------------------
# Median intensity as built-in predictor
# ---------------------------------------------------------------------------

class TestMedianIntensityPredictor:
    """Verify that median intensity is correctly used as a built-in predictor."""

    def _get_scores(self, tmp_path, columns_dict, variance_predictor_cols,
                    ion_var_values, ion_med_int_values=None):
        pids = _precursor_ids()
        ion_idx = _ion_index(pids)
        ml_file = _write_ml_info(tmp_path, pids, columns_dict)
        ion_var = _make_ion_variance(ion_idx, ion_var_values)
        ion_med_int = (_make_ion_median_intensity(ion_idx, ion_med_int_values)
                       if ion_med_int_values is not None else None)
        result = aq_vp.load_variance_predictor_scores(
            ml_file, ["S1", "S2"],
            variance_predictor_cols=variance_predictor_cols,
            ion_index=ion_idx, ion_variance=ion_var,
            ion_median_intensity=ion_med_int,
        )
        assert result is not None
        return np.array([result[ion] for ion in ion_idx])

    def test_intensity_alone_improves_prediction(self, tmp_path):
        """Adding median intensity as predictor should improve correlation
        with true variance when intensity is the dominant signal."""
        rng = np.random.RandomState(7)
        true_var = np.array([10.0 - 0.4 * i for i in range(N)])
        intensity = np.array([float(i) for i in range(N)])
        col_noisy = list(rng.randn(N))

        scores_no_int = self._get_scores(
            tmp_path, {"noisy": col_noisy}, ["noisy"],
            true_var, ion_med_int_values=None,
        )
        scores_with_int = self._get_scores(
            tmp_path, {"noisy": col_noisy}, ["noisy"],
            true_var, ion_med_int_values=intensity,
        )
        rho_no_int = abs(_spearman(true_var, scores_no_int))
        rho_with_int = abs(_spearman(true_var, scores_with_int))
        assert rho_with_int > rho_no_int + 0.1, (
            f"Expected intensity to improve prediction: "
            f"rho_with={rho_with_int:.3f}, rho_without={rho_no_int:.3f}"
        )

    def test_intensity_negative_coefficient(self, tmp_path):
        """Higher intensity → lower variance: model should assign negative
        weight to intensity and scores should track true variance."""
        intensity = np.array([float(i) for i in range(N)])
        true_var = np.array([float(N - 1 - i) * 0.5 for i in range(N)])
        col_vals = list(np.ones(N))  # uninformative quality col
        scores = self._get_scores(
            tmp_path, {"flat": col_vals}, ["flat"],
            true_var, ion_med_int_values=intensity,
        )
        rho = _spearman(true_var, scores)
        assert rho > 0.9

    def test_intensity_combined_with_quality_cols(self, tmp_path):
        """Intensity + opposing quality column: both should be recovered."""
        intensity = np.array([float(i) for i in range(N)])
        quality = np.array([float(N - 1 - i) for i in range(N)])
        true_var = 0.3 * intensity + 0.7 * quality
        scores = self._get_scores(
            tmp_path, {"quality": list(quality)}, ["quality"],
            true_var, ion_med_int_values=intensity,
        )
        rho = _spearman(true_var, scores)
        assert rho > 0.95

    def test_none_intensity_is_equivalent_to_no_intensity(self, tmp_path):
        """Passing None for ion_median_intensity should give the same results
        as not passing it at all."""
        col_vals = [float(i) for i in range(N)]
        true_var = [float(i) * 0.5 for i in range(N)]
        scores_none = self._get_scores(
            tmp_path, {"col": col_vals}, ["col"],
            true_var, ion_med_int_values=None,
        )
        pids = _precursor_ids()
        ion_idx = _ion_index(pids)
        ml_file = _write_ml_info(tmp_path, pids, {"col": col_vals})
        ion_var = _make_ion_variance(ion_idx, true_var)
        result = aq_vp.load_variance_predictor_scores(
            ml_file, ["S1", "S2"],
            variance_predictor_cols=["col"],
            ion_index=ion_idx, ion_variance=ion_var,
        )
        scores_default = np.array([result[ion] for ion in ion_idx])
        np.testing.assert_array_almost_equal(scores_none, scores_default)


# ---------------------------------------------------------------------------
# _fit_and_predict edge cases
# ---------------------------------------------------------------------------

class TestFitAndPredictEdgeCases:
    def test_returns_none_with_too_few_valid_ions(self):
        prec_features = pd.DataFrame({"col": [1.0, 2.0]}, index=["P1", "P2"])
        precursor_ids = np.array(["P1", "P2"])
        ion_index = pd.Index(["P1_FRGION_y1", "P2_FRGION_y1"])
        ion_var = pd.Series([0.1, 0.2], index=ion_index)
        result = aq_vp._fit_and_predict(
            prec_features, ["col"], precursor_ids, ion_index, ion_var
        )
        assert result is None

    def test_missing_features_get_median_fallback(self, tmp_path):
        """Ions whose precursor isn't in the ml_info_table should get the
        median predicted score instead of NaN."""
        pids = _precursor_ids()
        ion_idx = _ion_index(pids)
        col_vals = [float(i) for i in range(N)]
        true_var = [float(i) for i in range(N)]
        # Only write ml_info for half the precursors
        half_pids = pids[:N // 2]
        half_cols = col_vals[:N // 2]
        ml_file = _write_ml_info(tmp_path, half_pids, {"col": half_cols})
        ion_var = _make_ion_variance(ion_idx, true_var)
        result = aq_vp.load_variance_predictor_scores(
            ml_file, ["S1", "S2"], variance_predictor_cols=["col"],
            ion_index=ion_idx, ion_variance=ion_var,
        )
        assert result is not None
        # All scores should be finite (no NaN)
        assert all(np.isfinite(v) for v in result.values())


# ---------------------------------------------------------------------------
# Ion split pattern tests
# ---------------------------------------------------------------------------

class TestIonSplitPattern:
    def test_splits_frgion(self):
        parts = aq_vp._ION_SPLIT_PAT.split("PEP1_MOD1_2_FRGION_y3")
        assert parts[0] == "PEP1_MOD1_2"

    def test_splits_ms1isotopes(self):
        parts = aq_vp._ION_SPLIT_PAT.split("PEP1_MOD1_2_MS1ISOTOPES_0")
        assert parts[0] == "PEP1_MOD1_2"

    def test_no_split_for_precursor_ids(self):
        parts = aq_vp._ION_SPLIT_PAT.split("PEP1_MOD1_2")
        assert len(parts) == 1
        assert parts[0] == "PEP1_MOD1_2"
