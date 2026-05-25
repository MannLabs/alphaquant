import numpy as np
import pandas as pd
import pytest
import re

import alphaquant.diffquant.intensity_summarization as aq_summ
from alphaquant.cluster.cluster_ions import REGEX_FRGIONS_ISOTOPES


# ---------------------------------------------------------------------------
# Helpers — realistic ion names following the SEQ_..._MOD_..._CHARGE_... pattern
# ---------------------------------------------------------------------------

PROT = "PROT1"

# Precursor 1: charge 2
FRGION_Y3 = "SEQ_PEP1_MOD_MOD1_CHARGE_2_FRGION_y3_noloss_1"
FRGION_Y4 = "SEQ_PEP1_MOD_MOD1_CHARGE_2_FRGION_y4_noloss_1"
MS1ISO_0 = "SEQ_PEP1_MOD_MOD1_CHARGE_2_MS1ISOTOPES_0"
MS1ISO_1 = "SEQ_PEP1_MOD_MOD1_CHARGE_2_MS1ISOTOPES_1"
PREC_2 = "SEQ_PEP1_MOD_MOD1_CHARGE_2_PRECURSOR_2"

# Precursor 2: charge 3 (same peptide, different charge)
FRGION_C3_Y5 = "SEQ_PEP1_MOD_MOD1_CHARGE_3_FRGION_y5_noloss_1"

ALL_IONS = [FRGION_Y3, FRGION_Y4, MS1ISO_0, MS1ISO_1, PREC_2, FRGION_C3_Y5]
PEP2PROT = {ion: PROT for ion in ALL_IONS}


def _make_df(ions, values_per_sample):
    """Build a log2-intensity DataFrame from a dict {sample: [values_per_ion]}."""
    return pd.DataFrame(values_per_sample, index=ions)


# ---------------------------------------------------------------------------
# Tree building
# ---------------------------------------------------------------------------

class TestBuildTreeFromIonNames:
    def test_gene_root(self):
        tree = aq_summ.build_tree_from_ion_names(PROT, ALL_IONS)
        assert tree.name == PROT
        assert tree.type == "gene"

    def test_leaf_count(self):
        tree = aq_summ.build_tree_from_ion_names(PROT, ALL_IONS)
        assert len(tree.leaves) == len(ALL_IONS)

    def test_frgion_nodes_exist(self):
        tree = aq_summ.build_tree_from_ion_names(PROT, ALL_IONS)
        import anytree
        frgion_nodes = anytree.findall(tree, filter_=lambda n: n.type == "frgion")
        assert len(frgion_nodes) == 2  # one per charge state

    def test_ms1_nodes_exist(self):
        tree = aq_summ.build_tree_from_ion_names(PROT, ALL_IONS)
        import anytree
        ms1_nodes = anytree.findall(tree, filter_=lambda n: n.type == "ms1_isotopes")
        assert len(ms1_nodes) == 1

    def test_single_ion(self):
        tree = aq_summ.build_tree_from_ion_names(PROT, [FRGION_Y3])
        assert len(tree.leaves) == 1


# ---------------------------------------------------------------------------
# Grouping logic
# ---------------------------------------------------------------------------

class TestComputeSummarizationGroups:
    def test_empty_summarization_nodes(self):
        groups, remaining = aq_summ.compute_summarization_groups(PEP2PROT, ALL_IONS, [])
        assert groups == []
        assert remaining == set(ALL_IONS)

    def test_frgion_only(self):
        groups, remaining = aq_summ.compute_summarization_groups(PEP2PROT, ALL_IONS, ["frgion"])
        summarized_leaves = set()
        for _name, leaves, _prot in groups:
            summarized_leaves.update(leaves)
        # All FRGION base ions should be summarized
        assert FRGION_Y3 in summarized_leaves
        assert FRGION_Y4 in summarized_leaves
        assert FRGION_C3_Y5 in summarized_leaves
        # MS1 and precursor should remain
        assert MS1ISO_0 in remaining
        assert MS1ISO_1 in remaining
        assert PREC_2 in remaining
        # Two groups: one per charge state
        assert len(groups) == 2

    def test_ms1_only(self):
        groups, remaining = aq_summ.compute_summarization_groups(PEP2PROT, ALL_IONS, ["ms1_isotopes"])
        assert len(groups) == 1
        assert FRGION_Y3 in remaining
        assert FRGION_Y4 in remaining

    def test_frgion_and_ms1(self):
        groups, remaining = aq_summ.compute_summarization_groups(PEP2PROT, ALL_IONS, ["frgion", "ms1_isotopes"])
        summarized_leaves = set()
        for _name, leaves, _prot in groups:
            summarized_leaves.update(leaves)
        assert FRGION_Y3 in summarized_leaves
        assert MS1ISO_0 in summarized_leaves
        # Only precursor should remain
        assert remaining == {PREC_2}

    def test_mod_seq_charge_splits_by_ion_type(self):
        groups, remaining = aq_summ.compute_summarization_groups(PEP2PROT, ALL_IONS, ["mod_seq_charge"])
        # Two mod_seq_charge nodes (charge 2 and charge 3).
        # Charge 2 has frgion + ms1 + precursor -> 3 groups.
        # Charge 3 has only frgion -> 1 group.
        assert len(groups) == 4
        assert len(remaining) == 0


# ---------------------------------------------------------------------------
# Naming — summarized names must be parseable by downstream tree builder
# ---------------------------------------------------------------------------

class TestSummarizedNamesParseable:
    """Verify that summarized ion names are matched by the REGEX_FRGIONS_ISOTOPES
    patterns, so the downstream tree builder can incorporate them."""

    def _level0_matches(self, name):
        for pattern, _node_type in REGEX_FRGIONS_ISOTOPES[0]:
            if re.match(pattern, name):
                return True
        return False

    def test_frgion_sum_name_parseable(self):
        groups, _ = aq_summ.compute_summarization_groups(PEP2PROT, [FRGION_Y3, FRGION_Y4], ["frgion"])
        for name, _leaves, _prot in groups:
            assert self._level0_matches(name), f"'{name}' not parseable by level-0 regex"

    def test_ms1_sum_name_parseable(self):
        groups, _ = aq_summ.compute_summarization_groups(PEP2PROT, [MS1ISO_0, MS1ISO_1], ["ms1_isotopes"])
        for name, _leaves, _prot in groups:
            assert self._level0_matches(name), f"'{name}' not parseable by level-0 regex"

    def test_mod_seq_charge_sum_names_parseable(self):
        groups, _ = aq_summ.compute_summarization_groups(PEP2PROT, ALL_IONS, ["mod_seq_charge"])
        for name, _leaves, _prot in groups:
            assert self._level0_matches(name), f"'{name}' not parseable by level-0 regex"


# ---------------------------------------------------------------------------
# DataFrame summarization
# ---------------------------------------------------------------------------

class TestSummarizeConditionDf:
    @pytest.fixture
    def simple_df(self):
        return _make_df(
            [FRGION_Y3, FRGION_Y4, MS1ISO_0],
            {"s1": np.log2([100.0, 50.0, 500.0]), "s2": np.log2([200.0, 80.0, 600.0])},
        )

    def test_frgion_sum_values(self, simple_df):
        groups = [
            ("SEQ_PEP1_MOD_MOD1_CHARGE_2_FRGION_SUM", [FRGION_Y3, FRGION_Y4], PROT),
        ]
        remaining = {MS1ISO_0}
        result = aq_summ.summarize_condition_df(simple_df, groups, remaining)

        assert "SEQ_PEP1_MOD_MOD1_CHARGE_2_FRGION_SUM" in result.index
        assert MS1ISO_0 in result.index
        # sum(100+50)=150 for s1
        assert np.isclose(result.loc["SEQ_PEP1_MOD_MOD1_CHARGE_2_FRGION_SUM", "s1"], np.log2(150.0))
        # sum(200+80)=280 for s2
        assert np.isclose(result.loc["SEQ_PEP1_MOD_MOD1_CHARGE_2_FRGION_SUM", "s2"], np.log2(280.0))

    def test_ms1_unchanged(self, simple_df):
        groups = [
            ("SUMMED", [FRGION_Y3, FRGION_Y4], PROT),
        ]
        remaining = {MS1ISO_0}
        result = aq_summ.summarize_condition_df(simple_df, groups, remaining)

        assert np.isclose(result.loc[MS1ISO_0, "s1"], simple_df.loc[MS1ISO_0, "s1"])

    def test_all_nan_stays_nan(self):
        df = _make_df(
            [FRGION_Y3, FRGION_Y4],
            {"s1": [np.nan, np.nan], "s2": np.log2([100.0, 50.0])},
        )
        groups = [("SUM", [FRGION_Y3, FRGION_Y4], PROT)]
        result = aq_summ.summarize_condition_df(df, groups, set())

        assert np.isnan(result.loc["SUM", "s1"])
        assert np.isclose(result.loc["SUM", "s2"], np.log2(150.0))

    def test_partial_nan_sums_available(self):
        df = _make_df(
            [FRGION_Y3, FRGION_Y4],
            {"s1": [np.log2(100.0), np.nan], "s2": np.log2([100.0, 50.0])},
        )
        groups = [("SUM", [FRGION_Y3, FRGION_Y4], PROT)]
        result = aq_summ.summarize_condition_df(df, groups, set())

        # Only y3 contributes in s1
        assert np.isclose(result.loc["SUM", "s1"], np.log2(100.0))

    def test_empty_group_skipped(self):
        df = _make_df([MS1ISO_0], {"s1": [np.log2(500.0)]})
        groups = [("SUM", [FRGION_Y3, FRGION_Y4], PROT)]  # neither present in df
        remaining = {MS1ISO_0}
        result = aq_summ.summarize_condition_df(df, groups, remaining)

        assert len(result) == 1
        assert MS1ISO_0 in result.index


# ---------------------------------------------------------------------------
# End-to-end: apply_summarization
# ---------------------------------------------------------------------------

class TestApplySummarization:
    @pytest.fixture
    def condition_dfs(self):
        ions = [FRGION_Y3, FRGION_Y4, MS1ISO_0, MS1ISO_1, PREC_2]
        df_c1 = _make_df(ions, {
            "s1": np.log2([100.0, 50.0, 500.0, 200.0, 800.0]),
            "s2": np.log2([120.0, 60.0, 520.0, 210.0, 820.0]),
        })
        df_c2 = _make_df(ions, {
            "s3": np.log2([90.0, 40.0, 480.0, 190.0, 790.0]),
            "s4": np.log2([110.0, 55.0, 510.0, 205.0, 810.0]),
        })
        pep2prot = {ion: PROT for ion in ions}
        return df_c1, df_c2, pep2prot

    def test_no_summarization(self, condition_dfs):
        df_c1, df_c2, pep2prot = condition_dfs
        r1, r2, rp = aq_summ.apply_summarization(df_c1, df_c2, pep2prot, [])
        assert r1 is df_c1  # unchanged object
        assert r2 is df_c2

    def test_frgion_reduces_row_count(self, condition_dfs):
        df_c1, df_c2, pep2prot = condition_dfs
        r1, r2, rp = aq_summ.apply_summarization(df_c1, df_c2, pep2prot, ["frgion"])
        # 5 ions -> 1 frgion sum + 2 ms1 + 1 precursor = 4
        assert len(r1) == 4
        assert len(r2) == 4

    def test_pep2prot_updated(self, condition_dfs):
        df_c1, df_c2, pep2prot = condition_dfs
        _, _, rp = aq_summ.apply_summarization(df_c1, df_c2, pep2prot, ["frgion"])
        # Every row in the result should have a protein mapping
        for ion in set(rp.keys()):
            assert rp[ion] == PROT

    def test_frgion_and_ms1(self, condition_dfs):
        df_c1, df_c2, pep2prot = condition_dfs
        r1, r2, rp = aq_summ.apply_summarization(df_c1, df_c2, pep2prot, ["frgion", "ms1_isotopes"])
        # 1 frgion sum + 1 ms1 sum + 1 precursor = 3
        assert len(r1) == 3

    def test_asymmetric_conditions(self):
        """Ion present in c1 but not c2 — only ions in both conditions are summed."""
        df_c1 = _make_df([FRGION_Y3, FRGION_Y4], {"s1": np.log2([100.0, 50.0])})
        df_c2 = _make_df([FRGION_Y3], {"s2": np.log2([90.0])})
        pep2prot = {FRGION_Y3: PROT, FRGION_Y4: PROT}

        r1, r2, rp = aq_summ.apply_summarization(df_c1, df_c2, pep2prot, ["frgion"])
        # Y4 is absent from c2, so only Y3 (present in both) is used for the sum
        assert len(r1) == 1
        assert len(r2) == 1
        assert np.isclose(r1.iloc[0, 0], np.log2(100.0))
        assert np.isclose(r2.iloc[0, 0], np.log2(90.0))
