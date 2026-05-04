import alphaquant.cluster.cluster_utils as aq_clust_clusterutils
import anytree



def test_find_node_parent_at_level():

    udo = anytree.Node("Udo", type = 'granddad')
    marc = anytree.Node("Marc", parent=udo, type = 'dad')
    lian = anytree.Node("Lian", parent=marc, type = 'base')
    dan = anytree.Node("Dan", parent=udo, type ='dad')
    jet = anytree.Node("Jet", parent=dan, type ='base')
    jan = anytree.Node("Jan", parent=dan, type ='base')
    joe = anytree.Node("Joe", parent=dan, type ='base')

    assert aq_clust_clusterutils.find_node_parent_at_level(lian, 'dad').name == 'Marc'
    assert aq_clust_clusterutils.find_node_parent_at_level(lian, 'granddad').name == 'Udo'
    assert aq_clust_clusterutils.find_node_parent_at_level(jet, 'dad').name == 'Dan'

test_find_node_parent_at_level()


def test_check_if_node_is_included():

    udo = anytree.Node("Udo", type = 'gene', cluster =-1)
    marc = anytree.Node("Marc", parent=udo, type = 'dad', cluster = 0)
    lian = anytree.Node("Lian", parent=marc, type = 'base', cluster = 0)
    dan = anytree.Node("Dan", parent=udo, type ='dad', cluster = 1)
    jet = anytree.Node("Jet", parent=dan, type ='base', cluster = 0)
    jan = anytree.Node("Jan", parent=dan, type ='base', cluster = 0)
    joe = anytree.Node("Joe", parent=dan, type ='base', cluster = 0)
    assert aq_clust_clusterutils.check_if_node_is_included(jet) == False
    assert aq_clust_clusterutils.check_if_node_is_included(lian) == True
    assert aq_clust_clusterutils.check_if_node_is_included(dan) == False
    assert aq_clust_clusterutils.check_if_node_is_included(marc) == True




def test_remove_unnecessary_attributes():

    # Prepare a small tree with attributes
    root = anytree.Node("root", fcs="root_fcs", another_attr="root_other")
    child1 = anytree.Node("child1", parent=root, fcs="child1_fcs", another_attr="child1_other")
    child2 = anytree.Node("child2", parent=root, fcs="child2_fcs", another_attr="child2_other")
    # Leaf nodes
    leaf1 = anytree.Node("leaf1", parent=child1, fcs="leaf1_fcs")
    leaf2 = anytree.Node("leaf2", parent=child2, fcs="leaf2_fcs")

    # Invoke the spell to remove the 'fcs' attribute
    aq_clust_clusterutils.remove_unnecessary_attributes(root, ["fcs"])

    # Check if 'fcs' attribute is removed from all nodes
    all_nodes = [root, child1, child2, leaf1, leaf2]
    for node in all_nodes:
        assert not hasattr(node, "fcs"), f"Attribute 'fcs' was not removed from {node.name}."

    # Check if 'another_attr' remained untouched in nodes where it existed
    nodes_with_another_attr = [root, child1, child2]
    for node in nodes_with_another_attr:
        assert hasattr(node, "another_attr"), f"Attribute 'another_attr' was wrongly removed from {node.name}."

    print("All tests passed!")







def test_iterate_through_tree_levels_bottom_to_top():
    root = anytree.Node("root", level=0)
    child1 = anytree.Node("child1", parent=root, level=1)
    child2 = anytree.Node("child2", parent=root, level=1)
    child1_1 = anytree.Node("child1_1", parent=child1, level=2)
    child2_1 = anytree.Node("child2_1", parent=child2, level=2)

    expected_levels = [
        ["child1_1", "child2_1"],  # Level 2 nodes
        ["child1", "child2"],      # Level 1 nodes
        ["root"]                   # Level 0 node
    ]

    for index, nodes in enumerate(aq_clust_clusterutils.iterate_through_tree_levels_bottom_to_top(root)):
        actual_level_node_names = [node.name for node in nodes]
        assert actual_level_node_names == expected_levels[index], f"Level {index} does not match expected nodes."
        print(f"Level {index} matches expected nodes.")


def test_remove_outlier_fragion_childs_complete():
    """
    Comprehensive test that remove_outlier_fragion_childs() correctly filters fragments.
    Tests both normal mode and PTM mode.
    """
    import alphaquant.config.variables as aqvariables

    # Save original PTM_FRAGMENT_SELECTION setting
    original_ptm_setting = aqvariables.PTM_FRAGMENT_SELECTION

    try:
        # ========== Test Normal Mode (PTM_FRAGMENT_SELECTION = False) ==========
        aqvariables.PTM_FRAGMENT_SELECTION = False
        print("\n" + "="*60)
        print("TESTING NORMAL MODE (PTM_FRAGMENT_SELECTION = False)")
        print("="*60)

        # Test 1: With 4 or fewer fragments, all should be kept
        print("\n=== Test 1: 4 fragments (all kept) ===")
        fragments_4 = [
            anytree.Node("frag1", z_val=-2.0),
            anytree.Node("frag2", z_val=-0.5),
            anytree.Node("frag3", z_val=0.5),
            anytree.Node("frag4", z_val=2.0),
        ]

        result_4 = aq_clust_clusterutils.remove_outlier_fragion_childs(fragments_4)
        assert len(result_4) == 4, f"Expected 4 fragments, got {len(result_4)}"
        result_names = [f.name for f in result_4]
        expected_names = [f.name for f in fragments_4]
        assert set(result_names) == set(expected_names), "All 4 fragments should be kept"
        print(f"✓ All 4 fragments kept: {result_names}")

        # Verify is_outlier_fragment flags
        for frag in fragments_4:
            assert hasattr(frag, 'is_outlier_fragment'), "is_outlier_fragment flag should be set"
            assert frag.is_outlier_fragment == False, f"{frag.name} should not be marked as outlier"
        print(f"✓ is_outlier_fragment flags correctly set (all False)")

        # Test 2: With >4 fragments, only 4 middle ones (by z-value) should be kept
        # (median_idx - 2 : median_idx + 2 gives 4 elements due to Python slicing)
        print("\n=== Test 2: 7 fragments (4 middle kept) ===")
        fragments_7 = [
            anytree.Node("frag1", z_val=-3.0),  # Most negative (excluded)
            anytree.Node("frag2", z_val=-2.0),  # Should be kept
            anytree.Node("frag3", z_val=-1.0),  # Should be kept
            anytree.Node("frag4", z_val=0.0),   # Should be kept (median)
            anytree.Node("frag5", z_val=1.0),   # Should be kept
            anytree.Node("frag6", z_val=2.0),   # Excluded
            anytree.Node("frag7", z_val=3.0),   # Most positive (excluded)
        ]

        result_7 = aq_clust_clusterutils.remove_outlier_fragion_childs(fragments_7)
        assert len(result_7) == 4, f"Expected 4 fragments, got {len(result_7)}"

        result_zvals = sorted([f.z_val for f in result_7])
        expected_zvals = [-2.0, -1.0, 0.0, 1.0]
        assert result_zvals == expected_zvals, f"Expected z-vals {expected_zvals}, got {result_zvals}"
        print(f"✓ Correct 4 middle fragments kept with z-values: {result_zvals}")

        # Verify is_outlier_fragment flags
        assert fragments_7[0].is_outlier_fragment == True, "frag1 (z=-3.0) should be marked as outlier"
        assert fragments_7[5].is_outlier_fragment == True, "frag6 (z=2.0) should be marked as outlier"
        assert fragments_7[6].is_outlier_fragment == True, "frag7 (z=3.0) should be marked as outlier"
        for i in range(1, 5):
            assert fragments_7[i].is_outlier_fragment == False, f"frag{i+1} should not be marked as outlier"
        print(f"✓ is_outlier_fragment flags correctly set (3 outliers, 4 inliers)")

        # Test 3: With 10 fragments
        print("\n=== Test 3: 10 fragments (4 middle kept) ===")
        fragments_10 = [
            anytree.Node(f"frag{i}", z_val=float(i-5)) for i in range(10)
        ]
        # z_vals: -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0

        result_10 = aq_clust_clusterutils.remove_outlier_fragion_childs(fragments_10)
        assert len(result_10) == 4, f"Expected 4 fragments, got {len(result_10)}"

        result_zvals = sorted([f.z_val for f in result_10])
        expected_zvals = [-2.0, -1.0, 0.0, 1.0]
        assert result_zvals == expected_zvals, f"Expected z-vals {expected_zvals}, got {result_zvals}"
        print(f"✓ Correct 4 middle fragments kept with z-values: {result_zvals}")

        # Verify is_outlier_fragment flags
        expected_outliers = [0, 1, 2, 7, 8, 9]  # Indices of outliers
        expected_inliers = [3, 4, 5, 6]   # Indices of inliers
        for idx in expected_outliers:
            assert fragments_10[idx].is_outlier_fragment == True, f"frag{idx} should be marked as outlier"
        for idx in expected_inliers:
            assert fragments_10[idx].is_outlier_fragment == False, f"frag{idx} should not be marked as outlier"
        print(f"✓ is_outlier_fragment flags correctly set (6 outliers, 4 inliers)")

        # ========== PTM flag should not change this generic classic helper ==========
        aqvariables.PTM_FRAGMENT_SELECTION = True
        fragments_ptm_flag = [
            anytree.Node(f"frag{i}", z_val=float(i-10)) for i in range(20)
        ]
        result_ptm_flag = aq_clust_clusterutils.remove_outlier_fragion_childs(
            fragments_ptm_flag)
        assert len(result_ptm_flag) == 4

        print("\n=== All tests passed! ===")

    finally:
        # Restore original setting
        aqvariables.PTM_FRAGMENT_SELECTION = original_ptm_setting



# ---------------------------------------------------------------------------
# Tests for combine_zvalues / sum_and_re_scale_zvalues (new on this branch)
# ---------------------------------------------------------------------------
import numpy as np
import pytest


class TestCombineZvalues:
    """Tests for combine_zvalues dispatch and individual modes."""

    def test_single_value_returns_unchanged(self):
        assert aq_clust_clusterutils.combine_zvalues([1.5], mode="stouffer_decorrelation") == 1.5
        assert aq_clust_clusterutils.combine_zvalues([1.5], mode="mean_z") == 1.5
        assert aq_clust_clusterutils.combine_zvalues([1.5], mode="median_z") == 1.5
        assert aq_clust_clusterutils.combine_zvalues([1.5], mode="min_median_max_z") == 1.5

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError, match="Unknown aggregation mode"):
            aq_clust_clusterutils.combine_zvalues([1.0, 2.0], mode="bogus")

    def test_mean_z(self):
        zvals = [1.0, 2.0, 3.0]
        result = aq_clust_clusterutils.combine_zvalues(zvals, mode="mean_z")
        assert result == pytest.approx(2.0)

    def test_median_z(self):
        zvals = [1.0, 5.0, 2.0]
        result = aq_clust_clusterutils.combine_zvalues(zvals, mode="median_z")
        assert result == pytest.approx(2.0)

    def test_min_median_max_z_small_input(self):
        """n <= 3 should fall back to Stouffer on all values."""
        zvals = np.array([1.0, 2.0])
        result_mmm = aq_clust_clusterutils.combine_zvalues(zvals, mode="min_median_max_z")
        result_stouffer = aq_clust_clusterutils.sum_and_re_scale_zvalues(zvals, rho=0.0)
        assert result_mmm == pytest.approx(result_stouffer)

    def test_min_median_max_z_large_input(self):
        """n > 3 should use the 3-point summary."""
        zvals = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        result = aq_clust_clusterutils.combine_zvalues(zvals, mode="min_median_max_z")
        expected = aq_clust_clusterutils.sum_and_re_scale_zvalues(
            np.array([-2.0, 0.0, 2.0]), rho=0.0
        )
        assert result == pytest.approx(expected)

    def test_stouffer_decorrelation_mode_delegates(self):
        zvals = np.array([1.0, 2.0, 3.0])
        result = aq_clust_clusterutils.combine_zvalues(zvals, rho=0.0, mode="stouffer_decorrelation")
        expected = aq_clust_clusterutils.sum_and_re_scale_zvalues(zvals, rho=0.0)
        assert result == pytest.approx(expected)

    def test_stouffer_icc_remains_own_legacy_mode(self):
        zvals = np.array([1.0, 2.0, 3.0])
        result = aq_clust_clusterutils.combine_zvalues(zvals, rho=0.5, mode="stouffer_icc")
        expected = aq_clust_clusterutils.sum_and_re_scale_zvalues(zvals, rho=0.5)
        assert result == pytest.approx(expected)


class TestSumAndReScaleZvaluesWithRho:
    """Tests for the rho (ICC) parameter added to sum_and_re_scale_zvalues."""

    def test_rho_zero_is_classic_stouffer(self):
        zvals = np.array([1.0, 1.0, 1.0])
        result = aq_clust_clusterutils.sum_and_re_scale_zvalues(zvals, rho=0.0)
        assert result > 1.0  # combining concordant evidence should amplify

    def test_higher_rho_gives_more_conservative_result(self):
        zvals = np.array([2.0, 2.0, 2.0, 2.0])
        z_rho0 = aq_clust_clusterutils.sum_and_re_scale_zvalues(zvals, rho=0.0)
        z_rho05 = aq_clust_clusterutils.sum_and_re_scale_zvalues(zvals, rho=0.5)
        assert abs(z_rho05) < abs(z_rho0)

    def test_rho_one_returns_same_as_single_value(self):
        """With perfect correlation (rho=1), DEFF=n so sigma=n, reducing
        sum/sqrt(n*n) = mean, which for identical values is the value itself."""
        zvals = np.array([1.5, 1.5, 1.5])
        result = aq_clust_clusterutils.sum_and_re_scale_zvalues(zvals, rho=1.0)
        assert result == pytest.approx(1.5, abs=0.05)

    def test_single_value_ignores_rho(self):
        result = aq_clust_clusterutils.sum_and_re_scale_zvalues(np.array([2.5]), rho=0.8)
        assert result == pytest.approx(2.5)


def test_aggregate_node_properties_ignores_residual_decorrelation_exclusions():
    root = anytree.Node("parent", type="frgion", is_included=True, cluster=0)
    anytree.Node(
        "keep1",
        parent=root,
        type="base",
        is_included=True,
        cluster=0,
        z_val=1.0,
        fc=1.0,
        cv=0.2,
        min_intensity=10.0,
        total_intensity=10.0,
        min_reps=2,
        fraction_consistent=1.0,
    )
    anytree.Node(
        "drop",
        parent=root,
        type="base",
        is_included=True,
        cluster=0,
        z_val=100.0,
        fc=100.0,
        cv=0.2,
        min_intensity=10.0,
        total_intensity=10.0,
        min_reps=2,
        fraction_consistent=1.0,
        exclude_residual_decorrelation=True,
    )
    anytree.Node(
        "keep2",
        parent=root,
        type="base",
        is_included=True,
        cluster=0,
        z_val=3.0,
        fc=3.0,
        cv=0.4,
        min_intensity=20.0,
        total_intensity=20.0,
        min_reps=4,
        fraction_consistent=1.0,
    )

    aq_clust_clusterutils.aggregate_node_properties(root, only_use_mainclust=True, peptide_outlier_filtering=False)

    expected_z = aq_clust_clusterutils.combine_zvalues(np.array([1.0, 3.0]), rho=0.0, mode="stouffer_decorrelation")
    assert root.z_val == expected_z
    assert root.fc == 2.0
    assert root.min_reps == 3.0


def test_apply_ptm_fragment_selection_after_residual_decorrelation():
    root = anytree.Node("protein", type="gene")
    frgion = anytree.Node(
        "frgion",
        parent=root,
        type="frgion",
        is_included=True,
        cluster=0,
    )
    zvals = [10.0, -0.1, 3.0, 0.2, -7.0]
    children = []
    for idx, z_val in enumerate(zvals):
        children.append(anytree.Node(
            f"base{idx}",
            parent=frgion,
            type="base",
            is_included=True,
            cluster=0,
            z_val=z_val,
        ))
    children[0].exclude_residual_decorrelation = True

    dropped, parents = aq_clust_clusterutils.apply_ptm_fragment_selection(
        [root], max_keep=2)

    assert dropped == 2
    assert parents == 1
    assert children[0].exclude_residual_decorrelation is True
    assert not getattr(children[1], "exclude_ptm_fragment_selection", False)
    assert getattr(children[2], "exclude_ptm_fragment_selection", False)
    assert not getattr(children[3], "exclude_ptm_fragment_selection", False)
    assert getattr(children[4], "exclude_ptm_fragment_selection", False)


if __name__ == "__main__":
    test_remove_outlier_fragion_childs_complete()
    print("\n" + "="*60)
    print("ALL TESTS COMPLETED SUCCESSFULLY!")
    print("="*60)
