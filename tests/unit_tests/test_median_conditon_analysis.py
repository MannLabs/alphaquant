import alphaquant.multicond.median_condition_analysis as aq_multicond_mca
import anytree
import numpy as np
import pandas as pd
import pytest


""" test ProteoformConditionAligner """

def simulate_list_of_protein_nodes_from_cluster_matrix(cluster_matrix):
    list_of_protein_nodes = []
    num_condition = len(cluster_matrix[0])
    for cond_idx in range(num_condition):
        list_of_cluster_idxs_of_every_peptide = [x[cond_idx] for x in cluster_matrix]
        list_of_protein_nodes.append(simulate_protein_node_w_clustermatrix_subset(list_of_cluster_idxs_of_every_peptide, f"cond{cond_idx}"))
    return list_of_protein_nodes 
        

def simulate_protein_node_w_clustermatrix_subset(list_of_cluster_idxs_of_every_peptide, condition):
    condition = anytree.Node([condition, "median_reference"])
    protein = anytree.Node("protein", parent=condition)
    for idx, cluster in enumerate(list_of_cluster_idxs_of_every_peptide):
        peptide = anytree.Node(f"pep_{idx}", parent=protein)
        peptide.cluster = cluster
        peptide.fc = cluster
    return protein



def test_ProteoformConditionAligner1():
    cluster_matrix = [[0, 0, 1, 0],
                    [1, 1, 0, 1],
                    [0, 0, 1, 1],
                    [1, 1, 0, 1],
                    [1, 2, 3, 4]
                    ]
    expected_number_of_proteoforms = 3
    expected_peptide_groups = ["pep_0;pep_1;pep_3", "pep_2", "pep_4"]
    simulated_nodes = simulate_list_of_protein_nodes_from_cluster_matrix(cluster_matrix)
    aligner = aq_multicond_mca.ProteoformConditionAligner(simulated_nodes)
    assert len(set(aligner.proteoform_df["proteoform_id"])) == expected_number_of_proteoforms
    assert set(aligner.proteoform_df["peptides"]) == set(expected_peptide_groups)


def test_ProteoformConditionAligner2():
    cluster_matrix = [[0, 1, 1, 0], 
                    [1, 1, 0, 1], 
                    [0, 0, 1, 0], 
                    [1, 1, 0, 1]]
    expected_number_of_proteoforms = 2
    expected_peptide_groups = ["pep_0", "pep_1;pep_2;pep_3"]
    simulated_nodes = simulate_list_of_protein_nodes_from_cluster_matrix(cluster_matrix)
    aligner = aq_multicond_mca.ProteoformConditionAligner(simulated_nodes)
    assert len(set(aligner.proteoform_df["proteoform_id"])) == expected_number_of_proteoforms
    assert set(aligner.proteoform_df["peptides"]) == set(expected_peptide_groups)


""" test ProteoformPeptideDfCreator """


def create_proteoform_peptide_df_creator(list_of_fcs, list_of_numpep_per_cluster, number_conditions):
    list_of_conditions = [f"cond_{i}" for i in range(number_conditions)]
    nodes_same_protein_different_conditions = [
        create_protein_node(list_of_fcs, list_of_numpep_per_cluster, condition) 
        for condition in list_of_conditions
    ]
    return aq_multicond_mca.ProteoformPeptideDfCreator(nodes_same_protein_different_conditions)

def create_protein_node(list_of_fcs, list_of_numpep_per_cluster, condition):
    condition_node = anytree.Node([condition, "median_reference"])
    protein_node = anytree.Node("protein", parent=condition_node)
    for cluster_idx, numpep in enumerate(list_of_numpep_per_cluster):
        for pep_idx in range(numpep):
            peptide_node = anytree.Node(f"pep_{cluster_idx}_{pep_idx}", parent=protein_node)
            peptide_node.cluster = cluster_idx
            peptide_node.fc = list_of_fcs[cluster_idx]
    return protein_node

def get_expected_cluster_values(list_of_numpep_per_cluster):
    return [
        cluster_idx
        for cluster_idx, numpep in enumerate(list_of_numpep_per_cluster)
        for _ in range(numpep)
    ]

def get_expected_fc_values(list_of_numpeps_per_cluster, list_of_fcs):
    return [
        list_of_fcs[cluster_idx]
        for cluster_idx, numpep in enumerate(list_of_numpeps_per_cluster)
        for _ in range(numpep)
    ]

def assert_dataframe_dimensions(creator, list_of_numpep_per_cluster, num_conditions):
    num_peptides = sum(list_of_numpep_per_cluster)
    assert creator.peptide_cluster_df.shape == (num_peptides, num_conditions), "Incorrect dimensions for peptide_cluster_df"
    assert creator.peptide_fc_df.shape == (num_peptides, num_conditions), "Incorrect dimensions for peptide_fc_df"

def assert_cluster_values(creator, list_of_numpep_per_cluster):
    expected_cluster_values = get_expected_cluster_values(list_of_numpep_per_cluster)
    for col in creator.peptide_cluster_df.columns:
        assert all(creator.peptide_cluster_df[col] == expected_cluster_values), f"Mismatch in clusters for {col}"

def assert_fc_values(creator, list_of_numpeps_per_cluster, list_of_fcs):
    expected_fcs = get_expected_fc_values(list_of_numpeps_per_cluster, list_of_fcs)
    for col in creator.peptide_fc_df.columns:
        assert all(creator.peptide_fc_df[col] == expected_fcs), f"Mismatch in fc values for {col}"

@pytest.mark.parametrize("list_of_fcs, list_of_numpep_per_cluster, number_conditions", [
    ([1.1, -1.1], [2, 4], 5),
    ([1.1, -1.1, 10, 3], [2, 4, 12, 200], 60),
    ([1.1], [1], 2)
])
def test_proteoform_peptide_df_creator(list_of_fcs, list_of_numpep_per_cluster, number_conditions):
    creator = create_proteoform_peptide_df_creator(list_of_fcs, list_of_numpep_per_cluster, number_conditions)
    
    assert_dataframe_dimensions(creator, list_of_numpep_per_cluster, number_conditions)
    assert_cluster_values(creator, list_of_numpep_per_cluster)
    assert_fc_values(creator, list_of_numpep_per_cluster, list_of_fcs)


""" test ProteoformDfCreator """


def create_proteoform_df_creator(groups_of_peptide_clusters, peptide_fc_df, protein_name):
    return aq_multicond_mca.ProteoformDfCreator(groups_of_peptide_clusters, peptide_fc_df, protein_name)

def assert_initialization(creator):
    assert creator.proteoform_df is not None, "proteoform_df should not be None"

def assert_proteoform_grouping(creator, groups_of_peptide_clusters, protein_name):
    for idx, group in enumerate(groups_of_peptide_clusters):
        proteoform_id = f"{protein_name}_{idx}"
        group_df = creator.proteoform_df[creator.proteoform_df['proteoform_id'] == proteoform_id]
        assert all(group_df['peptides'] == ';'.join(group)), "Incorrect peptides in proteoform"
        assert all(group_df['protein'] == protein_name), "Incorrect protein name"

def assert_fold_change_data(creator, peptide_fc_df):
    for peptide in peptide_fc_df.index:
        for cond in peptide_fc_df.columns:
            assert np.isclose(creator.proteoform_df.loc[peptide, cond], peptide_fc_df.loc[peptide, cond]), f"Mismatch in fold change data for {peptide}, {cond}"

def assert_dataframe_structure(creator, peptide_fc_df, groups_of_peptide_clusters):
    expected_columns = ['peptides', 'protein', 'proteoform_id'] + list(peptide_fc_df.columns)
    assert set(creator.proteoform_df.columns) == set(expected_columns), "Incorrect dataframe columns"
    expected_rows = sum(len(group) for group in groups_of_peptide_clusters)
    assert len(creator.proteoform_df) == expected_rows, "Incorrect number of rows in dataframe"

def generate_peptide_fc_df(groups_of_peptide_clusters, number_conditions):
    np.random.seed(42)  # for reproducibility
    peptide2conds = {
        peptide: [np.random.random() for _ in range(number_conditions)]
        for group in groups_of_peptide_clusters
        for peptide in group
    }
    peptide2conds_df = pd.DataFrame(peptide2conds).T
    peptide2conds_df.columns = [f"cond_{idx}" for idx in range(number_conditions)]
    return peptide2conds_df

@pytest.mark.parametrize("groups_of_peptide_clusters, number_conditions", [
    ([["A", "B", "C"], ["D", "E"]], 3),
    ([["A", "B"], ["D", "E"], ["C"], ["F", "G", "H"]], 3)
])
def test_proteoform_df_creator(groups_of_peptide_clusters, number_conditions):
    peptide_fc_df = generate_peptide_fc_df(groups_of_peptide_clusters, number_conditions)
    protein_name = "protein"
    
    creator = create_proteoform_df_creator(groups_of_peptide_clusters, peptide_fc_df, protein_name)

    assert_initialization(creator)
    assert_proteoform_grouping(creator, groups_of_peptide_clusters, protein_name)
    assert_fold_change_data(creator, peptide_fc_df)
    assert_dataframe_structure(creator, peptide_fc_df, groups_of_peptide_clusters)