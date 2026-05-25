import uuid
import pandas as pd
import numpy as np
import alphaquant.diffquant.background_distributions as aq_diff_bg
import pytest
import time
from unittest.mock import patch, MagicMock


def generate_random_input(num_pep, sample2cond_df):
    pepnames = list(map(lambda _idx: str(uuid.uuid4()), range(num_pep)))  # gives uuid strings for each peptide
    randarrays = 10 + 1.5 * np.random.randn(num_pep, sample2cond_df.shape[0])
    df_intens = pd.DataFrame(randarrays, columns=sample2cond_df["sample"].tolist())
    df_intens.insert(0, "peptides", pepnames)
    df_intens = df_intens.set_index("peptides")
    return df_intens

@pytest.fixture
def sample2cond_df():
    return pd.DataFrame({'sample': ['A1', 'A2', 'A3', 'B1', 'B2', 'B3'], 'condition': ['A', 'A', 'A', 'B', 'B', 'B']})

@pytest.fixture
def fixed_input(sample2cond_df):
    return generate_random_input(1000, sample2cond_df)

@pytest.fixture
def background_distributions(fixed_input):
    """Create background distributions for testing caching"""
    condbg = aq_diff_bg.ConditionBackgrounds(fixed_input, {})
    # Collect unique BackGroundDistribution objects from ion2background
    seen = set()
    bg_list = []
    for bg in condbg.ion2background.values():
        bg_id = id(bg)
        if bg_id not in seen:
            seen.add(bg_id)
            bg_list.append(bg)
        if len(bg_list) >= 5:
            break
    return bg_list

def test_condition_backgrounds(fixed_input):
    condbg = aq_diff_bg.ConditionBackgrounds(fixed_input, {})
    assert condbg.ion2background.keys() == condbg.ion2nonNanvals.keys()


class TestBackgroundDistributionCaching:
    """Test suite for the improved background distribution caching functionality"""

    def test_background_cache_key_generation(self, background_distributions):
        """Test that cache keys are generated correctly and consistently"""
        bg1, bg2 = background_distributions[0], background_distributions[1]

        # Test that cache keys are consistent
        key1 = bg1.get_cache_key()
        key2 = bg1.get_cache_key()
        assert key1 == key2, "Cache keys should be consistent for the same object"

        # Test that different backgrounds have different keys
        key_bg1 = bg1.get_cache_key()
        key_bg2 = bg2.get_cache_key()
        assert key_bg1 != key_bg2, "Different backgrounds should have different cache keys"

        # Test cache key structure
        assert isinstance(key_bg1, tuple), "Cache key should be a tuple"
        assert len(key_bg1) == 6, "Cache key should have 6 elements"

        # Test that all elements are hashable
        hash(key_bg1)  # Should not raise an exception

    def test_get_background_cache_key_symmetry(self, background_distributions):
        """Test that the cache key function is symmetric (order-independent)"""
        bg1, bg2 = background_distributions[0], background_distributions[1]

        # Test symmetry: _get_background_cache_key(bg1, bg2) == _get_background_cache_key(bg2, bg1)
        key1 = aq_diff_bg._get_background_cache_key(bg1, bg2)
        key2 = aq_diff_bg._get_background_cache_key(bg2, bg1)

        assert key1 == key2, "Cache key should be symmetric regardless of argument order"

    def test_get_subtracted_bg_caching(self, background_distributions):
        """Test that get_subtracted_bg properly caches results"""
        bg1, bg2 = background_distributions[0], background_distributions[1]
        cache_dict = {}
        p2z = {}

        # First call should create and cache the result
        with patch('alphaquant.diffquant.background_distributions.SubtractedBackgrounds') as mock_constructor:
            mock_subtracted_bg = MagicMock()
            mock_constructor.return_value = mock_subtracted_bg

            result1 = aq_diff_bg.get_subtracted_bg(cache_dict, bg1, bg2, p2z)

            # Should have called the constructor once
            assert mock_constructor.call_count == 1
            assert result1 == mock_subtracted_bg

            # Cache should contain one entry
            assert len(cache_dict) == 1

        # Second call should use cached result
        with patch('alphaquant.diffquant.background_distributions.SubtractedBackgrounds') as mock_constructor:
            result2 = aq_diff_bg.get_subtracted_bg(cache_dict, bg1, bg2, p2z)

            # Should NOT have called the constructor
            assert mock_constructor.call_count == 0
            assert result2 == result1  # Should return the same cached object

            # Cache should still contain only one entry
            assert len(cache_dict) == 1

    def test_get_subtracted_bg_caching_symmetry(self, background_distributions):
        """Test that caching works symmetrically (same cache hit for swapped arguments)"""
        bg1, bg2 = background_distributions[0], background_distributions[1]
        cache_dict = {}
        p2z = {}

        # Create initial cache entry
        result1 = aq_diff_bg.get_subtracted_bg(cache_dict, bg1, bg2, p2z)

        # Test that swapped arguments hit the same cache
        with patch('alphaquant.diffquant.background_distributions.SubtractedBackgrounds') as mock_constructor:
            result2 = aq_diff_bg.get_subtracted_bg(cache_dict, bg2, bg1, p2z)

            # Should NOT create a new object (cache hit)
            assert mock_constructor.call_count == 0
            assert result2 == result1

            # Cache should still contain only one entry
            assert len(cache_dict) == 1

    def test_get_doublediff_bg_caching(self, background_distributions):
        """Test that get_doublediff_bg properly caches results"""
        bg1, bg2 = background_distributions[0], background_distributions[1]
        cache_dict = {}
        p2z = {}

        # First call should create and cache the result
        with patch('alphaquant.diffquant.background_distributions.SubtractedBackgrounds') as mock_constructor:
            mock_subtracted_bg = MagicMock()
            mock_constructor.return_value = mock_subtracted_bg

            result1 = aq_diff_bg.get_doublediff_bg(bg1, bg2, cache_dict, p2z)

            # Should have called the constructor once
            assert mock_constructor.call_count == 1
            assert result1 == mock_subtracted_bg

        # Second call should use cached result
        with patch('alphaquant.diffquant.background_distributions.SubtractedBackgrounds') as mock_constructor:
            result2 = aq_diff_bg.get_doublediff_bg(bg1, bg2, cache_dict, p2z)

            # Should NOT have called the constructor
            assert mock_constructor.call_count == 0
            assert result2 == result1

    def test_cache_performance_improvement(self, background_distributions):
        """Test that caching provides significant performance improvement"""
        bg1, bg2 = background_distributions[0], background_distributions[1]
        p2z = {}

        # Measure time for first call (no cache)
        cache_dict = {}
        start_time = time.time()
        result1 = aq_diff_bg.get_subtracted_bg(cache_dict, bg1, bg2, p2z)
        first_call_time = time.time() - start_time

        # Measure time for second call (with cache)
        start_time = time.time()
        result2 = aq_diff_bg.get_subtracted_bg(cache_dict, bg1, bg2, p2z)
        second_call_time = time.time() - start_time

        # Cached call should be significantly faster (at least 10x)
        assert second_call_time < first_call_time / 10, \
            f"Cached call ({second_call_time:.6f}s) should be much faster than first call ({first_call_time:.6f}s)"

        # Results should be identical
        assert result1 == result2

    def test_cache_with_many_different_backgrounds(self, background_distributions):
        """Test caching behavior with multiple different background pairs"""
        cache_dict = {}
        p2z = {}
        results = []

        # Create cache entries for all combinations
        for i in range(len(background_distributions)):
            for j in range(i + 1, len(background_distributions)):
                bg1, bg2 = background_distributions[i], background_distributions[j]
                result = aq_diff_bg.get_subtracted_bg(cache_dict, bg1, bg2, p2z)
                results.append(result)

        expected_cache_size = len(background_distributions) * (len(background_distributions) - 1) // 2
        assert len(cache_dict) == expected_cache_size, \
            f"Cache should contain {expected_cache_size} entries, but has {len(cache_dict)}"

        # Verify all results are different (different background pairs should produce different results)
        unique_results = set(id(result) for result in results)
        assert len(unique_results) == len(results), "All background pairs should produce unique results"

    def test_cache_key_uniqueness_across_different_distributions(self, fixed_input):
        """Test that cache keys are unique across different background distributions"""
        condbg1 = aq_diff_bg.ConditionBackgrounds(fixed_input, {})
        condbg2 = aq_diff_bg.ConditionBackgrounds(fixed_input, {})  # Different instance with same data

        # Get one background from each
        bg1_from_condbg1 = next(iter(condbg1.ion2background.values()))
        bg1_from_condbg2 = next(iter(condbg2.ion2background.values()))

        key1 = bg1_from_condbg1.get_cache_key()
        key2 = bg1_from_condbg2.get_cache_key()

        # They should be equal if they represent the same range and data
        assert key1 == key2, "Background distributions with identical parameters should have the same cache key"

    def test_cache_key_efficiency_and_reliability(self, background_distributions):
        """Test that cache keys provide reliable and consistent caching compared to string representations"""
        bg1, bg2 = background_distributions[0], background_distributions[1]

        # Test consistency - new cache keys should always be identical for same inputs
        keys = [aq_diff_bg._get_background_cache_key(bg1, bg2) for _ in range(100)]
        assert all(key == keys[0] for key in keys), "Cache keys should be consistent across multiple calls"

        # Test that our cache key is symmetric (major advantage over string approach)
        key1 = aq_diff_bg._get_background_cache_key(bg1, bg2)
        key2 = aq_diff_bg._get_background_cache_key(bg2, bg1)
        assert key1 == key2, "Cache key should be symmetric (order-independent)"

        # Test hashability and immutability (essential for dict keys)
        cache_key = aq_diff_bg._get_background_cache_key(bg1, bg2)
        assert isinstance(cache_key, tuple), "Cache key should be a tuple (immutable)"
        hash(cache_key)  # Should not raise an exception

        # Test that cache key uses stable properties
        # Each element should be a tuple containing numeric properties
        assert len(cache_key) == 2, "Cache key should have two tuple elements (one per background)"
        for bg_key in cache_key:
            assert isinstance(bg_key, tuple), "Each background key should be a tuple"
            assert len(bg_key) == 6, "Each background key should have 6 elements"


# ---------------------------------------------------------------------------
# Tests for _split_by_ion_type, _has_multiple_ion_types, and split background building
# (new on add_summarization_approach branch)
# ---------------------------------------------------------------------------
import alphaquant.config.variables as aqvariables


class TestSplitByIonType:
    def test_correct_masks(self):
        index = pd.Index([
            "PEP1_FRGION_y3",
            "PEP1_FRGION_y4",
            "PEP2_MS1ISOTOPES_0",
            "PEP2_MS1ISOTOPES_1",
        ])
        result = aq_diff_bg.ConditionBackgrounds._split_by_ion_type(index)
        assert list(result["FRGION"]) == [True, True, False, False]
        assert list(result["MS1ISOTOPES"]) == [False, False, True, True]

    def test_all_frgion(self):
        index = pd.Index(["A_FRGION_y1", "B_FRGION_b2"])
        result = aq_diff_bg.ConditionBackgrounds._split_by_ion_type(index)
        assert result["FRGION"].all()
        assert not result["MS1ISOTOPES"].any()

    def test_empty_index(self):
        index = pd.Index([])
        result = aq_diff_bg.ConditionBackgrounds._split_by_ion_type(index)
        assert len(result["FRGION"]) == 0
        assert len(result["MS1ISOTOPES"]) == 0


class TestHasMultipleIonTypes:
    def setup_method(self):
        self._original_config = aqvariables.CONFIG_DICT

    def teardown_method(self):
        aqvariables.CONFIG_DICT = self._original_config

    def test_true_when_both_present(self):
        aqvariables.CONFIG_DICT = {
            "ion_hierarchy": {"fragion": {}, "ms1iso": {}}
        }
        assert aq_diff_bg.ConditionBackgrounds._has_multiple_ion_types() is True

    def test_false_when_only_fragion(self):
        aqvariables.CONFIG_DICT = {
            "ion_hierarchy": {"fragion": {}}
        }
        assert aq_diff_bg.ConditionBackgrounds._has_multiple_ion_types() is False

    def test_false_when_config_is_none(self):
        aqvariables.CONFIG_DICT = None
        assert aq_diff_bg.ConditionBackgrounds._has_multiple_ion_types() is False


class TestConditionBackgroundsSplit:
    """Integration-level tests: verify split vs. single-pool mode."""

    def _make_mixed_df(self, n_frg=50, n_ms1=50, n_samples=4):
        rng = np.random.RandomState(42)
        frg_names = [f"PEP{i}_FRGION_y{j}" for i in range(n_frg) for j in [1]]
        ms1_names = [f"PEP{i}_MS1ISOTOPES_{j}" for i in range(n_ms1) for j in [0]]
        all_names = frg_names + ms1_names
        data = 10 + rng.randn(len(all_names), n_samples)
        cols = [f"S{i}" for i in range(n_samples)]
        return pd.DataFrame(data, index=all_names, columns=cols)

    def setup_method(self):
        self._original_config = aqvariables.CONFIG_DICT

    def teardown_method(self):
        aqvariables.CONFIG_DICT = self._original_config

    def test_split_mode_assigns_all_ions(self):
        aqvariables.CONFIG_DICT = {
            "ion_hierarchy": {"fragion": {}, "ms1iso": {}}
        }
        df = self._make_mixed_df()
        cb = aq_diff_bg.ConditionBackgrounds(df, {}, split_by_ion_type=True)
        assert set(cb.ion2background.keys()) == set(df.index)

    def test_single_pool_mode_assigns_all_ions(self):
        df = self._make_mixed_df()
        cb = aq_diff_bg.ConditionBackgrounds(df, {}, split_by_ion_type=False)
        assert set(cb.ion2background.keys()) == set(df.index)

    def test_split_false_when_config_missing(self):
        aqvariables.CONFIG_DICT = None
        df = self._make_mixed_df()
        cb = aq_diff_bg.ConditionBackgrounds(df, {}, split_by_ion_type=True)
        assert set(cb.ion2background.keys()) == set(df.index)
