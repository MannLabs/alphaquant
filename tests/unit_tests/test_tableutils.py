import pandas as pd
import alphaquant.tables.tableutils as aq_table_tableutils
import numpy as np

def test_init_with_ml_score():
    sample_df_ml_score = pd.DataFrame({
        'ml_score': [0.8, 0.3, 0.5, 0.1, 0.9]
    })
    normalizer = aq_table_tableutils.QualityScoreNormalizer(sample_df_ml_score)
    assert 'quality_score' in normalizer.results_df.columns
    assert 'ml_score' not in normalizer.results_df.columns
    print("test_init_with_ml_score passed")

def test_init_with_consistency_score():
    sample_df_consistency = pd.DataFrame({
        'consistency_score': [70, 85, 60, 95, 80]
    })
    normalizer = aq_table_tableutils.QualityScoreNormalizer(sample_df_consistency)
    assert 'quality_score' in normalizer.results_df.columns
    assert 'consistency_score' not in normalizer.results_df.columns
    print("test_init_with_consistency_score passed")

def test_init_with_invalid_df():
    invalid_df = pd.DataFrame({'invalid_column': [1, 2, 3]})
    try:
        aq_table_tableutils.QualityScoreNormalizer(invalid_df)
    except ValueError:
        print("test_init_with_invalid_df passed")
    else:
        assert False, "ValueError not raised"

def test_normalize_quality_score_ml():
    sample_df_ml_score = pd.DataFrame({
        'ml_score': [0.8, 0.3, 0.5, 0.1, 0.9]
    })
    normalizer = aq_table_tableutils.QualityScoreNormalizer(sample_df_ml_score)
    quality_scores = normalizer.results_df['quality_score'].values
    assert np.isclose(quality_scores.min(), 0.2)
    assert np.isclose(quality_scores.max(), 1.0)
    print("test_normalize_quality_score_ml passed")

def test_normalize_quality_score_consistency():
    sample_df_consistency = pd.DataFrame({
        'consistency_score': [70, 85, 60, 95, 80]
    })
    normalizer = aq_table_tableutils.QualityScoreNormalizer(sample_df_consistency)
    quality_scores = normalizer.results_df['quality_score'].values
    assert np.isclose(quality_scores.min(), 0.2)
    assert np.isclose(quality_scores.max(), 1.0)
    print("test_normalize_quality_score_consistency passed")

def test_rank_normalization_ordered_scores():
    # Test case: Ordered scores
    scores1 = pd.Series([1, 2, 3, 4, 5])
    normalized_ranks_descending = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores1, higher_is_better=True)
    print(normalized_ranks_descending)
    assert np.allclose(normalized_ranks_descending, [0.2, 0.4, 0.6, 0.8, 1.0])
    normalized_ranks_ascending = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores1, higher_is_better=False)
    assert np.allclose(normalized_ranks_ascending, [1.0, 0.8, 0.6, 0.4, 0.2])

def test_rank_normalization_unordered_scores():
    scores2 = pd.Series([3, 1, 4, 5, 2])
    normalized_ranks_descending_unordered = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores2, higher_is_better=True)
    assert np.allclose(normalized_ranks_descending_unordered, [0.6, 0.2, 0.8, 1.0, 0.4])
    normalized_ranks_ascending_unordered = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores2, higher_is_better=False)
    assert np.allclose(normalized_ranks_ascending_unordered, [0.6, 1.0, 0.4, 0.2, 0.8])

def test_rank_normalization_negative_scores():
    scores3 = pd.Series([-5, -3, -1, 1, 3])
    normalized_ranks_descending_negative = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores3, higher_is_better=True)
    assert np.allclose(normalized_ranks_descending_negative, [0.2, 0.4, 0.6, 0.8, 1.0])
    normalized_ranks_ascending_negative = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores3, higher_is_better=False)
    assert np.allclose(normalized_ranks_ascending_negative, [1.0, 0.8, 0.6, 0.4, 0.2])

def test_rank_normalization_identical_elements():
    scores4 = pd.Series([1, 1, 1, 1, 1])
    normalized_ranks_identical_higher = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores4, higher_is_better=True)
    assert np.allclose(normalized_ranks_identical_higher, [0.6, 0.6, 0.6, 0.6, 0.6])
    normalized_ranks_identical_lower = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores4, higher_is_better=False)
    assert np.allclose(normalized_ranks_identical_lower, [0.6, 0.6, 0.6, 0.6, 0.6])

def test_rank_normalization_tied_scores():
    scores5 = pd.Series([1, 2, 2, 3, 3])
    normalized_ranks_ties_higher = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores5, higher_is_better=True)
    print(normalized_ranks_ties_higher)
    assert np.allclose(normalized_ranks_ties_higher, [0.2, 0.5, 0.5, 0.9, 0.9])
    normalized_ranks_ties_lower = aq_table_tableutils.QualityScoreNormalizer._perform_rank_normalization(scores5, higher_is_better=False)
    assert np.allclose(normalized_ranks_ties_lower, [1.0, 0.7, 0.7, 0.3, 0.3])

    print("test_perform_rank_normalization passed")
