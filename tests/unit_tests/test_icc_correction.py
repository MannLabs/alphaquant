import numpy as np
import anytree
import pytest

import alphaquant.cluster.icc_correction as aq_icc


# ---------------------------------------------------------------------------
# Helpers for building synthetic trees
# ---------------------------------------------------------------------------

def _make_base_node(name, parent, z_val):
    return anytree.Node(name, parent=parent, type="base", z_val=z_val,
                        is_included=True, cluster=0)


def _make_protein_tree(gene_name, group_zvals, node_type="frgion", p_val=0.5):
    """Build a minimal protein tree.

    Args:
        gene_name: Name for the root (gene) node.
        group_zvals: list of lists — each inner list holds z-values for
            one group node's base children.
        node_type: Type string for the group-level nodes (e.g. "frgion").
        p_val: Gene-level p-value attached to the root.

    Returns:
        anytree.Node: Root of the protein tree.
    """
    root = anytree.Node(gene_name, type="gene", p_val=p_val,
                        is_included=True, cluster=-1)
    for i, zvals in enumerate(group_zvals):
        grp = anytree.Node(f"{gene_name}_grp{i}", parent=root,
                           type=node_type, is_included=True, cluster=0)
        for j, z in enumerate(zvals):
            _make_base_node(f"{gene_name}_grp{i}_base{j}", parent=grp, z_val=z)
    return root


# ---------------------------------------------------------------------------
# _compute_icc_from_tree
# ---------------------------------------------------------------------------

class TestComputeIccFromTree:
    def test_returns_none_when_no_matching_nodes(self):
        root = _make_protein_tree("P1", [[0.1, 0.2], [0.3, 0.4]], node_type="frgion")
        assert aq_icc._compute_icc_from_tree(root, "ms1_isotopes") is None

    def test_returns_none_when_too_few_groups(self):
        root = _make_protein_tree("P1", [[0.1, 0.2, 0.3]], node_type="frgion")
        assert aq_icc._compute_icc_from_tree(root, "frgion") is None

    def test_returns_none_when_too_few_ions(self):
        root = _make_protein_tree("P1", [[0.1], [0.2], [0.3]], node_type="frgion")
        assert aq_icc._compute_icc_from_tree(root, "frgion") is None

    def test_identical_within_group_gives_high_icc(self):
        group_zvals = [[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]]
        root = _make_protein_tree("P1", group_zvals, node_type="frgion")
        icc = aq_icc._compute_icc_from_tree(root, "frgion")
        assert icc is not None
        assert icc > 0.9

    def test_no_between_group_variance_gives_zero_icc(self):
        rng = np.random.RandomState(42)
        group_zvals = [list(rng.randn(5)) for _ in range(5)]
        # Shift all groups to the same mean
        for gv in group_zvals:
            m = np.mean(gv)
            for i in range(len(gv)):
                gv[i] -= m
        root = _make_protein_tree("P1", group_zvals, node_type="frgion")
        icc = aq_icc._compute_icc_from_tree(root, "frgion")
        assert icc is not None
        assert icc < 0.1

    def test_icc_between_zero_and_one(self):
        rng = np.random.RandomState(99)
        group_zvals = [list(rng.randn(4) + i) for i in range(4)]
        root = _make_protein_tree("P1", group_zvals, node_type="frgion")
        icc = aq_icc._compute_icc_from_tree(root, "frgion")
        assert icc is not None
        assert 0.0 <= icc <= 1.0


# ---------------------------------------------------------------------------
# _estimate_null_icc_distribution
# ---------------------------------------------------------------------------

class TestEstimateNullIccDistribution:
    def test_only_null_proteins_are_used(self):
        """Proteins with p_val <= threshold should be excluded."""
        sig = _make_protein_tree("sig", [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                                 p_val=0.01)
        null = _make_protein_tree("null", [[0.1, 0.1, 0.1], [0.2, 0.2, 0.2],
                                           [0.3, 0.3, 0.3]], p_val=0.5)
        null_iccs, median = aq_icc._estimate_null_icc_distribution(
            [sig, null], "frgion"
        )
        assert len(null_iccs) == 1

    def test_empty_when_no_null_proteins(self):
        sig = _make_protein_tree("sig", [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                                 p_val=0.01)
        null_iccs, median = aq_icc._estimate_null_icc_distribution(
            [sig], "frgion"
        )
        assert null_iccs == []
        assert median == 0.0


# ---------------------------------------------------------------------------
# _compute_and_assign_icc_for_protein
# ---------------------------------------------------------------------------

class TestComputeAndAssignIccForProtein:
    def test_assigns_icc_to_target_nodes(self):
        root = _make_protein_tree("P1", [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                                  node_type="frgion")
        result = aq_icc._compute_and_assign_icc_for_protein(root, "frgion", 0.5)
        assert result is not None
        for node in anytree.search.findall(root, filter_=lambda n: n.type == "frgion"):
            assert hasattr(node, "icc_correction")
            assert hasattr(node, "icc_is_fallback")

    def test_caps_at_null_median(self):
        root = _make_protein_tree("P1", [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                                  node_type="frgion")
        icc_median = 0.1
        result = aq_icc._compute_and_assign_icc_for_protein(root, "frgion", icc_median)
        assert result <= icc_median

    def test_fallback_when_insufficient_data(self):
        root = _make_protein_tree("P1", [[1.0]], node_type="frgion")
        icc_median = 0.3
        result = aq_icc._compute_and_assign_icc_for_protein(root, "frgion", icc_median)
        assert result == icc_median
        frgion_nodes = anytree.search.findall(root, filter_=lambda n: n.type == "frgion")
        for node in frgion_nodes:
            assert node.icc_is_fallback is True

    def test_returns_none_when_no_target_nodes(self):
        root = _make_protein_tree("P1", [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                                  node_type="frgion")
        result = aq_icc._compute_and_assign_icc_for_protein(root, "ms1_isotopes", 0.5)
        assert result is None


# ---------------------------------------------------------------------------
# _has_node_type
# ---------------------------------------------------------------------------

class TestHasNodeType:
    def test_true_when_present(self):
        root = _make_protein_tree("P1", [[0.1, 0.2]], node_type="frgion")
        assert aq_icc._has_node_type([root], "frgion") is True

    def test_false_when_absent(self):
        root = _make_protein_tree("P1", [[0.1, 0.2]], node_type="frgion")
        assert aq_icc._has_node_type([root], "ms1_isotopes") is False
