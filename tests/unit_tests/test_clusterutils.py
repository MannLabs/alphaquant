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





def test_traverse_and_add_included_leaves_anytree():
    # Constructing the tree
    root = anytree.Node("root",  is_included=True, cluster=0)
    node1 = anytree.Node("node1", parent=root, is_included=True, cluster=0)
    node2 = anytree.Node("node2", parent=root, is_included=True, cluster=0)
    leaf1 = anytree.Node("leaf1", parent=node1, is_included=True, cluster=0)
    leaf2 = anytree.Node("leaf2", parent=node1, is_included=False, cluster=1)
    leaf3 = anytree.Node("leaf3", parent=node2, is_included=True, cluster=0)

    list_of_included_leaves = []
    aq_clust_clusterutils.traverse_and_add_included_leaves(root, list_of_included_leaves)
    print(list_of_included_leaves)
    # Assert conditions
    assert leaf1 in list_of_included_leaves, "leaf1 is missing from the result."
    assert leaf3 in list_of_included_leaves, "leaf3 is missing from the result."
    assert len(list_of_included_leaves) == 2, "The number of included leaves is incorrect."


    root = anytree.Node("root",  is_included=True, cluster=0)
    node1 = anytree.Node("node1", parent=root, is_included=False, cluster=1)
    node2 = anytree.Node("node2", parent=root, is_included=True, cluster=0)
    leaf1 = anytree.Node("leaf1", parent=node1, is_included=True, cluster=0)
    leaf2 = anytree.Node("leaf2", parent=node1, is_included=False, cluster=1)
    leaf3 = anytree.Node("leaf3", parent=node2, is_included=True, cluster=0)

    list_of_included_leaves = []
    aq_clust_clusterutils.traverse_and_add_included_leaves(root, list_of_included_leaves)
    print(list_of_included_leaves)
    # Assert conditions
    assert leaf1  not in list_of_included_leaves, "leaf1 should be excluded"
    assert leaf3 in list_of_included_leaves, "leaf3 is missing from the result."
    assert len(list_of_included_leaves) == 1, "The number of included leaves is incorrect."

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



