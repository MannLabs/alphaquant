import anytree
import pytest

import alphaquant.cluster.cluster_missingval as aq_missingval


def _reset_global():
    """Reset the module-level test-level global before each test."""
    aq_missingval.MISSINGVAL_TEST_LEVEL = None


# ---------------------------------------------------------------------------
# Helper: build minimal tree structures for each scenario
# ---------------------------------------------------------------------------

def _tree_with_mod_seq_charge():
    """gene -> seq -> mod_seq -> mod_seq_charge -> frgion -> base"""
    root = anytree.Node("gene1", type="gene")
    seq = anytree.Node("seq1", parent=root, type="seq")
    mod = anytree.Node("mod1", parent=seq, type="mod_seq")
    msc = anytree.Node("msc1", parent=mod, type="mod_seq_charge")
    frg = anytree.Node("frg1", parent=msc, type="frgion")
    anytree.Node("base1", parent=frg, type="base")
    anytree.Node("base2", parent=frg, type="base")
    return root


def _tree_mod_seq_above_leaves():
    """gene -> seq -> mod_seq -> base (no mod_seq_charge)"""
    root = anytree.Node("gene1", type="gene")
    seq = anytree.Node("seq1", parent=root, type="seq")
    mod = anytree.Node("mod1", parent=seq, type="mod_seq")
    anytree.Node("base1", parent=mod, type="base")
    anytree.Node("base2", parent=mod, type="base")
    return root


def _tree_seq_above_leaves():
    """gene -> seq -> base (precursor-only, no mod_seq)"""
    root = anytree.Node("gene1", type="gene")
    seq = anytree.Node("seq1", parent=root, type="seq")
    anytree.Node("base1", parent=seq, type="base")
    anytree.Node("base2", parent=seq, type="base")
    return root


def _tree_gene_above_leaves():
    """gene -> base (simplest hierarchy)"""
    root = anytree.Node("gene1", type="gene")
    anytree.Node("base1", parent=root, type="base")
    anytree.Node("base2", parent=root, type="base")
    return root


# ---------------------------------------------------------------------------
# Tests for determine_missingval_test_level
# ---------------------------------------------------------------------------

class TestDetermineMissingvalTestLevel:
    def setup_method(self):
        _reset_global()

    def test_mod_seq_charge_tree(self):
        root = _tree_with_mod_seq_charge()
        aq_missingval.determine_missingval_test_level(root)
        assert aq_missingval.MISSINGVAL_TEST_LEVEL == "mod_seq_charge"

    def test_mod_seq_above_leaves(self):
        root = _tree_mod_seq_above_leaves()
        aq_missingval.determine_missingval_test_level(root)
        assert aq_missingval.MISSINGVAL_TEST_LEVEL == "base"

    def test_seq_above_leaves(self):
        root = _tree_seq_above_leaves()
        aq_missingval.determine_missingval_test_level(root)
        assert aq_missingval.MISSINGVAL_TEST_LEVEL == "base"

    def test_gene_above_leaves(self):
        root = _tree_gene_above_leaves()
        aq_missingval.determine_missingval_test_level(root)
        assert aq_missingval.MISSINGVAL_TEST_LEVEL == "base"

    def test_unexpected_structure_raises(self):
        root = anytree.Node("root", type="unknown_top")
        anytree.Node("leaf", parent=root, type="base")
        with pytest.raises(ValueError, match="Unexpected tree structure"):
            aq_missingval.determine_missingval_test_level(root)


class TestGetNodesToTest:
    """Tests for MissingValProtNodeCreator._get_nodes_to_test."""

    def setup_method(self):
        _reset_global()

    def test_returns_mod_seq_charge_nodes_when_present(self):
        root = _tree_with_mod_seq_charge()
        nodes = aq_missingval.MissingValProtNodeCreator._get_nodes_to_test(root)
        assert all(n.type == "mod_seq_charge" for n in nodes)

    def test_returns_leaves_for_mod_seq_tree(self):
        root = _tree_mod_seq_above_leaves()
        nodes = aq_missingval.MissingValProtNodeCreator._get_nodes_to_test(root)
        assert all(n.type == "base" for n in nodes)
        assert set(n.name for n in nodes) == {"base1", "base2"}

    def test_returns_leaves_for_gene_tree(self):
        root = _tree_gene_above_leaves()
        nodes = aq_missingval.MissingValProtNodeCreator._get_nodes_to_test(root)
        assert all(n.type == "base" for n in nodes)
