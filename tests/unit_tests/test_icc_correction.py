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


def _make_protein_tree(gene_name, group_zvals, node_type="frgion", p_val=0.5,
                       group_p_vals=None):
    """Build a minimal protein tree.

    Args:
        gene_name: Name for the root (gene) node.
        group_zvals: list of lists — each inner list holds z-values for
            one group node's base children.
        node_type: Type string for the group-level nodes (e.g. "frgion").
        p_val: Gene-level p-value attached to the root.
        group_p_vals: Optional list of p-values for the group nodes.
            If None, each group gets p_val=0.5.

    Returns:
        anytree.Node: Root of the protein tree.
    """
    root = anytree.Node(gene_name, type="gene", p_val=p_val,
                        is_included=True, cluster=-1)
    for i, zvals in enumerate(group_zvals):
        grp_p = group_p_vals[i] if group_p_vals is not None else 0.5
        grp = anytree.Node(f"{gene_name}_grp{i}", parent=root,
                           type=node_type, is_included=True, cluster=0,
                           p_val=grp_p)
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
        null_iccs, perm_iccs, median = aq_icc._estimate_null_icc_distribution(
            [sig, null], "frgion"
        )
        assert len(null_iccs) == 1
        assert len(perm_iccs) > 0

    def test_empty_when_no_null_proteins(self):
        sig = _make_protein_tree("sig", [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                                 p_val=0.01)
        null_iccs, perm_iccs, median = aq_icc._estimate_null_icc_distribution(
            [sig], "frgion"
        )
        assert null_iccs == []
        assert perm_iccs == []
        assert median == 0.0

    def test_permutation_null_near_zero(self):
        """Permutation ICC should be much lower than observed ICC for correlated data."""
        rng = np.random.RandomState(7)
        group_zvals = [list(rng.randn(5) + offset) for offset in range(5)]
        root = _make_protein_tree("null", group_zvals, p_val=0.5)
        null_iccs, perm_iccs, median = aq_icc._estimate_null_icc_distribution(
            [root], "frgion"
        )
        assert len(perm_iccs) > 0
        assert np.median(perm_iccs) < median

    def test_significant_group_nodes_excluded_from_null(self):
        """Group nodes with p_val <= threshold should be dropped during null estimation."""
        group_zvals = [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]]
        group_p_vals = [0.5, 0.5, 0.01, 0.5]
        root = _make_protein_tree("null", group_zvals, p_val=0.5,
                                  group_p_vals=group_p_vals)
        icc_filtered = aq_icc._compute_icc_from_tree(
            root, "frgion", node_p_val_threshold=0.1
        )
        icc_unfiltered = aq_icc._compute_icc_from_tree(root, "frgion")
        assert icc_filtered is not None
        assert icc_unfiltered is not None
        assert icc_filtered != icc_unfiltered

    def test_all_group_nodes_significant_returns_none(self):
        """If every group node is significant, the protein contributes nothing."""
        group_zvals = [[1, 1, 1], [2, 2, 2], [3, 3, 3]]
        group_p_vals = [0.01, 0.01, 0.01]
        root = _make_protein_tree("null", group_zvals, p_val=0.5,
                                  group_p_vals=group_p_vals)
        icc = aq_icc._compute_icc_from_tree(
            root, "frgion", node_p_val_threshold=0.1
        )
        assert icc is None


# ---------------------------------------------------------------------------
# _assign_icc_to_all_proteins
# ---------------------------------------------------------------------------

class TestAssignIccToAllProteins:
    def test_assigns_uniform_icc_to_all_nodes(self):
        root = _make_protein_tree("P1", [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                                  node_type="frgion")
        icc_median = 0.25
        n = aq_icc._assign_icc_to_all_proteins([root], "frgion", icc_median)
        assert n == 3
        for node in anytree.search.findall(root, filter_=lambda n: n.type == "frgion"):
            assert node.icc_correction == icc_median

    def test_returns_zero_when_no_matching_nodes(self):
        root = _make_protein_tree("P1", [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                                  node_type="frgion")
        n = aq_icc._assign_icc_to_all_proteins([root], "ms1_isotopes", 0.3)
        assert n == 0


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
