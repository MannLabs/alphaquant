from anytree import Node
import alphaquant.cluster.proteoform_statistics as aq_clust_pformpvals
from unittest.mock import patch
import numpy as np

def test_get_leavenames_from_node():
    # Create a tree structure
    root = Node("root")
    child1 = Node("child1", parent=root)
    child2 = Node("child2", parent=root)
    leaf1 = Node("leaf1", parent=child1)
    leaf2 = Node("leaf2", parent=child1)
    leaf3 = Node("leaf3", parent=child2)

    # Test with take_median_ion=False
    result = aq_clust_pformpvals._get_leavenames_from_node(child1, take_median_ion=False)
    assert set(result) == {"leaf1", "leaf2"}


    with patch('alphaquant.cluster.cluster_utils.select_middle_leafs', return_value=[leaf2]):
        result = aq_clust_pformpvals._get_leavenames_from_node(child1, take_median_ion=True)
        assert result == ["leaf2"]


def test_get_cluster2nodes():
    # Create a tree structure
    root = Node("root")
    child1 = Node("child1", parent=root)
    child2 = Node("child2", parent=root)
    grandchild1 = Node("grandchild1", parent=child1)
    grandchild2 = Node("grandchild2", parent=child2)
    grandchild3 = Node("grandchild3", parent=child2)

    # Create childnode2clust_ordered dictionary
    childnode2clust_ordered = {
        child1: 0,
        child2: 1,
        grandchild1: 0,
        grandchild2: 1,
        grandchild3: 2
    }

    # Call the function
    result = aq_clust_pformpvals._get_cluster2nodes(childnode2clust_ordered)

    # Assert the results
    assert len(result) == 3
    assert set(result[0]) == {child1, grandchild1}
    assert set(result[1]) == {child2, grandchild2}
    assert set(result[2]) == {grandchild3}

def test_get_cluster2ions():
    # Create a tree structure
    root = Node("root")
    child1 = Node("child1", parent=root, fc = 0.5)
    child2 = Node("child2", parent=root, fc = 1)
    leaf1 = Node("ion1", parent=child1)
    leaf2 = Node("ion2", parent=child2)
    leaf3 = Node("ion3", parent=child2)

    # Create cluster2nodes dictionary
    cluster2nodes = {
        0: [child1],
        1: [child2]
    }

    # Test with take_median_ions=False
    result = aq_clust_pformpvals._get_cluster2ions(cluster2nodes, take_median_ions=False, take_median_node=True)

    # Assert the results
    assert len(result) == 2
    assert result[0] == ["ion1"]
    assert set(result[1]) == {"ion2", "ion3"}

    # Test with take_median_ions=True
    with patch('alphaquant.cluster.cluster_utils.select_middle_leafs', side_effect=lambda leaves: [leaves[0]]):
        result = aq_clust_pformpvals._get_cluster2ions(cluster2nodes, take_median_ions=True, take_median_node=True)

        # Assert the results
        assert len(result) == 2
        assert result[0] == ["ion1"]
        assert result[1] == ["ion2"]  # Only one ion selected as median


def test_assign_proteoform_statisitcs():
    # Create a tree structure
    root = Node("root", type="seq")
    child1 = Node("child1", parent=root, type="seq", fc = 0.5)
    child2 = Node("child2", parent=root, type="seq", fc = 1)
    leaf1 = Node("ion1", parent=child1)
    leaf2 = Node("ion2", parent=child2)
    leaf3 = Node("ion3", parent=child2)

    # Create childnode2clust_ordered dictionary
    childnode2clust_ordered = {
        child1: 0,
        child2: 1
    }

    # Mock objects
    normed_c1, normed_c2 = object(), object()
    ion2diffDist, p2z, deedpair2doublediffdist = object(), object(), object()

    # Mock calc_doublediff_score function
    with patch('alphaquant.diffquant.doublediff_analysis.calc_doublediff_score', return_value=(1.5, 0.05)) as mock_calc:
        aq_clust_pformpvals.add_proteoform_statistics_to_nodes(
            childnode2clust_ordered, False, normed_c1, normed_c2,
            ion2diffDist, p2z, deedpair2doublediffdist
        )

        # Assert calc_doublediff_score was called with correct arguments
        mock_calc.assert_called_once_with(
            ions1=['ion1'], ions2=['ion2', 'ion3'],
            normed_c1=normed_c1, normed_c2=normed_c2,
            ion2diffDist=ion2diffDist, p2z=p2z,
            deedpair2doublediffdist=deedpair2doublediffdist
        )

    # Assert proteoform stats were annotated correctly
    assert np.isnan(child1.proteoform_fcfc)
    assert np.isnan(child1.proteoform_pval)
    assert child2.proteoform_fcfc == -1.5
    assert child2.proteoform_pval == 0.05