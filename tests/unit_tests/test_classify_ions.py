import pytest
from anytree import Node, AnyNode
import pandas as pd
import numpy as np
import alphaquant.classify.classify_ions_stacked as aq_classify_ions_stacked
import alphaquant.classify.classify_ions as aqclassify
from sklearn.ensemble import RandomForestRegressor

# Mock class for testing
class MLInputTableCreatorTest(aq_classify_ions_stacked.MLInputTableCreator):
    def __init__(self, precursors, acquisition_info_df, replace_nans=False, numeric_threshold=0.2):
        self._precursors = precursors
        self._acquisition_info_df = acquisition_info_df
        self._replace_nans = replace_nans
        self._numeric_threshold = numeric_threshold

        self._merged_df = None
        
        self.X = None
        self.y = None
        self.featurenames = None
        self.ionnames = None

        self._define_merged_df()
        self._define_ionnames()
        self._remove_non_numeric_columns_from_merged_df()
        self._define_featurenames()
        self._define_X()
        self._define_y()

    def _define_merged_df(self):
        data = {
            "quant_id": ["Q1", "Q2", "Q3", "Q4", "Q5", "Q6"],
            "feature1": [1.0, 2.0, 2, 3.0, 3, np.nan],
            "feature2": [4, 5.0, 6.0, 2, 7.0, np.nan],
            "feature3": [7.0, 8.0, 9.0, 10.0, 11.0, np.nan],
            "feature4": [np.nan, 13.0, 14.0, 15.0,  16.0, np.nan],
            "feature5": [17.0, 2, 18.0, 19.0, 20.0, np.nan]
        }
        self._merged_df = pd.DataFrame(data)

def create_precursors():
    gene_parent = AnyNode(name="GeneParent", level="gene", fc=0.5)
    precursors = [
        AnyNode(name=f"Q{i}", fc=i+0.5, parent=gene_parent, level="mod_seq_charge")
        for i in range(1, 7)
    ]
    return precursors

@pytest.fixture
def setup_precursor_selector():
    protein1 = Node("Protein1", level="gene", fc=0.5)
    protein2 = Node("Protein2", level="gene", fc=0.8)
    protein3 = Node("Protein3", level="gene", fc=0.8)
    
    precursors = [
        Node(f"Precursor{i}", parent=protein1 if i <= 2 else protein2 if i <= 5 else protein3, level="mod_seq_charge")
        for i in range(1, 8)
    ]
    
    selector = aq_classify_ions_stacked.PrecursorForTrainingSelector(
        protein_nodes=[protein1, protein2], min_num_precursors=3, prot_fc_cutoff=0.75
    )
    return selector, precursors

def test_precursor_selector(setup_precursor_selector):
    selector, precursors = setup_precursor_selector
    
    assert precursors[2] in selector.precursors_suitable_for_training
    assert precursors[3] in selector.precursors_suitable_for_training
    assert precursors[4] in selector.precursors_suitable_for_training
    assert precursors[0] in selector.precursors_not_suitable_for_training
    assert precursors[1] in selector.precursors_not_suitable_for_training

@pytest.fixture
def setup_ml_input_table_creator():
    precursors = create_precursors()
    acquisition_info_df = pd.DataFrame({"quant_id": [f"Q{i}" for i in range(1, 7)], "extra_info": [100, 200, 100, 200, 100, 200]})
    ml_creator = MLInputTableCreatorTest(precursors, acquisition_info_df, replace_nans=True)
    return ml_creator, precursors

def test_ml_input_table_creator(setup_ml_input_table_creator):
    ml_creator, precursors = setup_ml_input_table_creator
    
    assert ml_creator._merged_df is not None
    assert "feature1" in ml_creator._merged_df.columns
    assert len(ml_creator.y) == 6
    assert not np.isnan(ml_creator.X).any()
    assert len(ml_creator.ionnames) == 6
    assert len(ml_creator.ionnames) == len(ml_creator.y)
    assert len(ml_creator.featurenames) == ml_creator.X.shape[1]
    
    name2fc = {node.name: node.fc - node.parent.fc for node in precursors}
    assert np.allclose(ml_creator.y, [name2fc[ionname] for ionname in ml_creator.ionnames])

def test_iterative_cross_predict():
    np.random.seed(42)  # for reproducibility
    X = np.random.rand(200, 3)
    y = np.random.rand(200)
    ionnames = np.random.rand(200)
    select_idxs = np.random.randint(low=0, high=199, size=15)
    control_ionnames = ionnames[select_idxs].copy()
    
    y_test_all, _, ionnames_all, _ = aqclassify.random_forest_iterative_cross_predict(
        X, y, ionnames, 5, RandomForestRegressor()
    )
    
    idxs_ionnames = [x for x in range(len(ionnames_all)) if ionnames_all[x] in control_ionnames]
    control_idxs_y = [x for x in select_idxs if ionnames[x] in ionnames_all]
    control_ys = y[control_idxs_y]
    
    y_test_all_control = np.array(y_test_all)[idxs_ionnames]
    
    assert set(control_ys) == set(y_test_all_control)