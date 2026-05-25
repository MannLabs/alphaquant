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
    """Build a minimal protein tree with one grouping level.

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


def _make_deep_tree(gene_name, p_val=0.5, seed=42):
    """Build a multi-level tree matching the real hierarchy.

    gene → seq → mod_seq → mod_seq_charge → {frgion, ms1_isotopes} → base
    """
    rng = np.random.RandomState(seed)
    root = anytree.Node(gene_name, type="gene", p_val=p_val,
                        is_included=True, cluster=-1)

    for s in range(3):
        seq = anytree.Node(f"{gene_name}_seq{s}", parent=root,
                           type="seq", is_included=True, cluster=0,
                           p_val=0.5, z_val=rng.randn())
        for m in range(2):
            mod = anytree.Node(f"{gene_name}_seq{s}_mod{m}", parent=seq,
                               type="mod_seq", is_included=True, cluster=0,
                               p_val=0.5, z_val=rng.randn())
            for c in range(2):
                msc = anytree.Node(f"{gene_name}_seq{s}_mod{m}_ch{c}",
                                   parent=mod, type="mod_seq_charge",
                                   is_included=True, cluster=0,
                                   p_val=0.5, z_val=rng.randn())
                frg = anytree.Node(f"{gene_name}_seq{s}_mod{m}_ch{c}_frg",
                                   parent=msc, type="frgion",
                                   is_included=True, cluster=0,
                                   p_val=0.5, z_val=rng.randn())
                for b in range(3):
                    _make_base_node(
                        f"{gene_name}_seq{s}_mod{m}_ch{c}_frg_b{b}",
                        parent=frg, z_val=rng.randn())

                ms1 = anytree.Node(f"{gene_name}_seq{s}_mod{m}_ch{c}_ms1",
                                   parent=msc, type="ms1_isotopes",
                                   is_included=True, cluster=0,
                                   p_val=0.5, z_val=rng.randn())
                for b in range(2):
                    _make_base_node(
                        f"{gene_name}_seq{s}_mod{m}_ch{c}_ms1_b{b}",
                        parent=ms1, z_val=rng.randn())
    return root


# ---------------------------------------------------------------------------
# _compute_icc_from_tree  (works at any level now)
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

    def test_works_with_non_base_children(self):
        """ICC at mod_seq_charge level uses ion-type children (not base)."""
        root = _make_deep_tree("P1")
        icc = aq_icc._compute_icc_from_tree(root, "mod_seq_charge")
        assert icc is None or 0.0 <= icc <= 1.0


# ---------------------------------------------------------------------------
# _collect_group_zvals  (generalized for any level)
# ---------------------------------------------------------------------------

class TestCollectGroupZvals:
    def test_collects_base_children_for_frgion(self):
        group_zvals = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
        root = _make_protein_tree("P1", group_zvals, node_type="frgion")
        result = aq_icc._collect_group_zvals(root, "frgion")
        assert len(result) == 3
        np.testing.assert_array_equal(result[0], [1.0, 2.0])

    def test_collects_non_base_children_for_higher_levels(self):
        """At mod_seq_charge level, children are frgion/ms1_isotopes nodes."""
        root = _make_deep_tree("P1")
        result = aq_icc._collect_group_zvals(root, "mod_seq_charge")
        assert len(result) > 0
        for arr in result:
            assert len(arr) >= 1

    def test_gene_level_returns_empty_due_to_single_group(self):
        """A single protein has 1 gene node → 1 group < _MIN_GROUPS → empty.

        Gene-level ICC uses the separate _estimate_gene_level_icc path.
        """
        root = _make_deep_tree("P1")
        result = aq_icc._collect_group_zvals(root, "gene")
        assert result == []

    def test_filters_by_p_val_threshold(self):
        group_zvals = [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]]
        group_p_vals = [0.5, 0.5, 0.01, 0.5]
        root = _make_protein_tree("P1", group_zvals, p_val=0.5,
                                  group_p_vals=group_p_vals)
        result = aq_icc._collect_group_zvals(root, "frgion",
                                             node_p_val_threshold=0.1)
        assert len(result) == 3

    def test_empty_when_all_filtered(self):
        group_zvals = [[1, 1, 1], [2, 2, 2], [3, 3, 3]]
        group_p_vals = [0.01, 0.01, 0.01]
        root = _make_protein_tree("P1", group_zvals, p_val=0.5,
                                  group_p_vals=group_p_vals)
        result = aq_icc._collect_group_zvals(root, "frgion",
                                             node_p_val_threshold=0.1)
        assert result == []


# ---------------------------------------------------------------------------
# _estimate_null_icc_distribution
# ---------------------------------------------------------------------------

class TestEstimateNullIccDistribution:
    def test_only_null_proteins_are_used(self):
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
        group_zvals = [[1.0, 1.1, 0.9], [2.0, 2.3, 1.7],
                       [3.0, 3.5, 2.5], [4.0, 3.8, 4.2]]
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
        group_zvals = [[1, 1, 1], [2, 2, 2], [3, 3, 3]]
        group_p_vals = [0.01, 0.01, 0.01]
        root = _make_protein_tree("null", group_zvals, p_val=0.5,
                                  group_p_vals=group_p_vals)
        icc = aq_icc._compute_icc_from_tree(
            root, "frgion", node_p_val_threshold=0.1
        )
        assert icc is None


# ---------------------------------------------------------------------------
# _estimate_gene_level_icc
# ---------------------------------------------------------------------------

class TestEstimateGeneLevelIcc:
    def _make_null_proteins(self, n_proteins=20, n_peptides_range=(3, 8),
                            p_val=0.5, seed=42):
        """Create a list of protein trees with seq children carrying z-values."""
        rng = np.random.RandomState(seed)
        protnodes = []
        for i in range(n_proteins):
            n_pep = rng.randint(n_peptides_range[0], n_peptides_range[1] + 1)
            root = anytree.Node(f"prot{i}", type="gene", p_val=p_val,
                                is_included=True, cluster=-1)
            for j in range(n_pep):
                anytree.Node(f"prot{i}_seq{j}", parent=root,
                             type="seq", is_included=True, cluster=0,
                             z_val=rng.randn())
            protnodes.append(root)
        return protnodes

    def test_returns_icc_distribution(self):
        protnodes = self._make_null_proteins(n_proteins=20)
        null_iccs, perm_iccs, median = aq_icc._estimate_gene_level_icc(protnodes)
        assert len(null_iccs) > 0
        assert len(perm_iccs) > 0
        assert 0.0 <= median <= 1.0

    def test_fallback_when_too_few_proteins(self):
        protnodes = self._make_null_proteins(n_proteins=3)
        null_iccs, perm_iccs, median = aq_icc._estimate_gene_level_icc(protnodes)
        assert null_iccs == []
        assert median == 0.0

    def test_excludes_significant_proteins(self):
        null_prots = self._make_null_proteins(n_proteins=20, p_val=0.5)
        sig_prots = self._make_null_proteins(n_proteins=5, p_val=0.01, seed=99)
        all_prots = null_prots + sig_prots
        null_iccs, _, median = aq_icc._estimate_gene_level_icc(all_prots)
        assert len(null_iccs) > 0

    def test_excludes_proteins_with_single_peptide(self):
        """Proteins with only 1 seq child should be skipped."""
        protnodes = []
        for i in range(20):
            root = anytree.Node(f"prot{i}", type="gene", p_val=0.5,
                                is_included=True, cluster=-1)
            anytree.Node(f"prot{i}_seq0", parent=root, type="seq",
                         is_included=True, cluster=0, z_val=0.1)
            protnodes.append(root)
        null_iccs, _, median = aq_icc._estimate_gene_level_icc(protnodes)
        assert null_iccs == []
        assert median == 0.0

    def test_shuffled_lower_than_observed_for_correlated_data(self):
        """When peptides within proteins are correlated, permutation ICC < observed."""
        rng = np.random.RandomState(77)
        protnodes = []
        for i in range(30):
            root = anytree.Node(f"prot{i}", type="gene", p_val=0.5,
                                is_included=True, cluster=-1)
            protein_effect = rng.randn() * 2
            for j in range(5):
                anytree.Node(f"prot{i}_seq{j}", parent=root, type="seq",
                             is_included=True, cluster=0,
                             z_val=protein_effect + rng.randn() * 0.3)
            protnodes.append(root)
        null_iccs, perm_iccs, median = aq_icc._estimate_gene_level_icc(protnodes)
        assert np.median(perm_iccs) < median


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

    def test_assigns_to_gene_node(self):
        root = _make_protein_tree("P1", [[1, 2], [3, 4], [5, 6]],
                                  node_type="frgion")
        n = aq_icc._assign_icc_to_all_proteins([root], "gene", 0.15)
        assert n == 1
        assert root.icc_correction == 0.15

    def test_assigns_to_higher_level_nodes(self):
        root = _make_deep_tree("P1")
        n = aq_icc._assign_icc_to_all_proteins([root], "mod_seq_charge", 0.1)
        msc_nodes = anytree.search.findall(root, filter_=lambda n: n.type == "mod_seq_charge")
        assert n == len(msc_nodes)
        for node in msc_nodes:
            assert node.icc_correction == 0.1


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

    def test_gene_type_found_at_root(self):
        root = _make_protein_tree("P1", [[0.1, 0.2]], node_type="frgion")
        assert aq_icc._has_node_type([root], "gene") is True

    def test_deep_tree_has_all_types(self):
        root = _make_deep_tree("P1")
        for node_type in ("gene", "seq", "mod_seq", "mod_seq_charge",
                          "frgion", "ms1_isotopes", "base"):
            assert aq_icc._has_node_type([root], node_type) is True
