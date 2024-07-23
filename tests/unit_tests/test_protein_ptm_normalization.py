import alphaquant.ptm.protein_ptm_normalization as aqptmnorm
import pandas as pd
import os
import shutil
import random
import numpy as np
import itertools
import glob
import pytest

random.seed(42)
np.random.seed(42)

RESULTS_DIR_SIMULATED_PROTEOME = "../../test_data/unit_tests/protein_ptm_normalization/toy_data/results_proteome"
RESULTS_DIR_SIMULATED_PTM = "../../test_data/unit_tests/protein_ptm_normalization/toy_data/results_ptm"

def test_p_value_damping():
    assert_that_half_fc_is_half_logp_value()
    assert_that_half_fc_is_half_logp_value_minus()
    assert_that_unsignificant_protein_retains_ptm_p_value()
    assert_that_no_protfc_is_ptmp_value()
    assert_that_same_fc_equals_p_value1()


def assert_that_half_fc_is_half_logp_value():
    reginfos = aqptmnorm.RegulationInfos(log2fc_ptm=1.5, p_value_ptm=1e-10, log2fc_protein=0.75,p_value_protein=0.001)
    p_valuedamper = aqptmnorm.PvalDamper(reginfos)
    damped_p_value = p_valuedamper.get_p_value()
    assert damped_p_value == 1e-5

def assert_that_half_fc_is_half_logp_value_minus():
    reginfos = aqptmnorm.RegulationInfos(log2fc_ptm=-1.5, p_value_ptm=1e-10, log2fc_protein=-0.75,p_value_protein=0.001)
    p_valuedamper = aqptmnorm.PvalDamper(reginfos)
    damped_p_value = p_valuedamper.get_p_value()
    assert damped_p_value == 1e-5

def assert_that_unsignificant_protein_retains_ptm_p_value():
    reginfos = aqptmnorm.RegulationInfos(log2fc_ptm=1.5, p_value_ptm=1e-10, log2fc_protein=0.75,p_value_protein=0.1)
    p_valuedamper = aqptmnorm.PvalDamper(reginfos)
    damped_p_value = p_valuedamper.get_p_value()
    assert damped_p_value == 1e-10

def assert_that_no_protfc_is_ptmp_value():
    reginfos =aqptmnorm.RegulationInfos(log2fc_ptm=1.5, p_value_ptm=1e-10, log2fc_protein=0,p_value_protein=0.001)
    p_valuedamper = aqptmnorm.PvalDamper(reginfos)
    damped_p_value = p_valuedamper.get_p_value()
    assert damped_p_value == 1e-10

def assert_that_same_fc_equals_p_value1():
    reginfos = aqptmnorm.RegulationInfos(log2fc_ptm=1.5, p_value_ptm=1e-10, log2fc_protein=1.5,p_value_protein=0.001)
    p_valuedamper = aqptmnorm.PvalDamper(reginfos)
    damped_p_value = p_valuedamper.get_p_value()
    assert damped_p_value == 1


"""Simulate example tables for testing the normalization"""


class ResultsDirSimulator():
    def __init__(self, results_dir_protein, results_dir_ptm, resultswriterconfig_vec):
        self._results_dir_protein = results_dir_protein
        self._results_dir_ptm = results_dir_ptm
        self._resultswriterconfig_vec = resultswriterconfig_vec
        self._create_toy_results_dirs()
        self._write_out_simulated_dataframes()

    
    def _create_toy_results_dirs(self):
        self._create_or_replace_folder(self._results_dir_protein)
        self._create_or_replace_folder(self._results_dir_ptm)
    
    @staticmethod
    def _create_or_replace_folder(folder):
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.makedirs(folder)

    def _write_out_simulated_dataframes(self):
        for idx, resultswriterconfig in enumerate(self._resultswriterconfig_vec):
            condpairname = f"condpair{idx}" #make up a name for the tables
            resultswriterconfig.filenames_w_config.append(condpairname)
            resutstable_simulator = ResultsTableSimulator(resultswriterconfig)
            self._save_dataframe(resutstable_simulator.protein_df, self._results_dir_protein, condpairname)
            self._save_dataframe(resutstable_simulator.ptm_df, self._results_dir_ptm, condpairname)

    @staticmethod
    def _save_dataframe(df, results_dir, name):
        df.to_csv(f"{results_dir}/{name}.results.tsv", sep = "\t", index = None)


class ResultsTableSimulator():
    def __init__(self, resultswriterconfig):
        self.protein_df = None
        self.ptm_df = None
        self._resultswriterconfig = resultswriterconfig
        self._simulate_dataframes_as_specified()
    
    def _simulate_dataframes_as_specified(self):
        protein_df_initial = self._simulate_protein_df()
        ptm_df_initial = self._simulate_ptm_df(protein_df_initial, offset_fc_to_protein=self._resultswriterconfig.offset_fc)
        self.protein_df =self._remove_rows_from_df(protein_df_initial, nrows_to_remove=self._resultswriterconfig.number_missing_proteins)
        self.ptm_df = self._remove_rows_from_df(ptm_df_initial, nrows_to_remove=self._resultswriterconfig.number_missing_ptms)

    def _simulate_protein_df(self):
        proteins = [f'name{x}' for x in range(self._resultswriterconfig.length)]
        fcs = self._get_fcs()
        p_value = [10**(-np.random.uniform(0, 10)) for x in range(self._resultswriterconfig.length)]
        return pd.DataFrame(data = {'protein' : proteins, 'log2fc' : fcs, 'p_value' : p_value})

    def _get_fcs(self):
        if self._resultswriterconfig.fc_for_all_proteins is not None:
            return [self._resultswriterconfig.fc_for_all_proteins for x in range( self._resultswriterconfig.length)]
        else:
            return [np.random.uniform(low = -2, high= 2) for x in range(self._resultswriterconfig.length)]
    
    @staticmethod
    def _simulate_ptm_df(protein_df, offset_fc_to_protein):
        ptm_df = protein_df.copy()
        ptm_df["swissprot"] = protein_df["protein"]
        ptm_df["protein"] = [f"{x}_protname_ptmsite" for x in protein_df["protein"]]
        ptm_df["log2fc"] = [x + offset_fc_to_protein for x in ptm_df["log2fc"]]
        return ptm_df
    
    @staticmethod
    def _remove_rows_from_df(results_df, nrows_to_remove):
        all_idxs = list(results_df.index)
        randomly_sampled_idxs_to_remove = random.sample(all_idxs, nrows_to_remove)
        return results_df.drop(index = randomly_sampled_idxs_to_remove)


class ResultsWriterConfig():
    def __init__(self, length, offset_fc, number_missing_proteins, number_missing_ptms, fc_for_all_proteins = None):
        self.length = length
        self.offset_fc = offset_fc
        self.number_missing_proteins = number_missing_proteins
        self.number_missing_ptms = number_missing_ptms
        self.filenames_w_config = []
        self.fc_for_all_proteins = fc_for_all_proteins



"""Check that normalized outputs are as expected from simulation"""





class PTMNormalizationTester:
    def __init__(self, simulate_new=True):
        self.results_config_vec = self.create_results_config_vec()
        self.filename2results_config = {}
        self.results_ptm = RESULTS_DIR_SIMULATED_PTM
        self.results_proteome = RESULTS_DIR_SIMULATED_PROTEOME
        self.simulate_results_dir(simulate_new)
        self.localizer = self.get_ptmtable_localizer()

    def create_results_config_vec(self):
        config1 = ResultsWriterConfig(50, 1, 8, 9)
        config2 = ResultsWriterConfig(50, -0.5, 8, 9)
        config3 = ResultsWriterConfig(50, +1.5, 38, 38)
        config4 = ResultsWriterConfig(50, 1.5, 0, 0, 0.5)
        results_config_vec = [config1] * 10 + [config2] * 10 + [config3] * 10 + [config4] * 2
        return results_config_vec

    def simulate_results_dir(self, simulate_new):
        if simulate_new:
            ResultsDirSimulator(self.results_proteome, self.results_ptm, self.results_config_vec)

    def get_ptmtable_localizer(self):
        return aqptmnorm.PTMtableLocalizer(self.results_ptm, self.results_proteome)

    def fcs_match_resultsconfig(self, log2fcs_list):
        for resconf in set(self.results_config_vec):
            comparevec = [resconf.offset_fc] * len(log2fcs_list)
            if np.allclose(log2fcs_list, comparevec):
                return True
        return False

@pytest.fixture(scope="module")
def ptm_normalization_tester():
    return PTMNormalizationTester()

def test_correct_localization(ptm_normalization_tester):
    ptmfile2name = ptm_normalization_tester.localizer.get_ptmfile2name()
    assert len(ptmfile2name) == len(ptm_normalization_tester.results_config_vec)
    
    list_of_filelists = [x.filenames_w_config for x in ptm_normalization_tester.results_config_vec]
    all_files = list(set(itertools.chain.from_iterable(list_of_filelists)))
    assert sorted(all_files) == sorted(list(ptmfile2name.values()))

def test_correct_normalization(ptm_normalization_tester):
    ptmnorm = aqptmnorm.PTMResultsNormalizer(ptm_normalization_tester.results_ptm, ptm_normalization_tester.results_proteome)
    results_dir_protnormed = ptmnorm.results_dir_protnormed
    
    all_files = glob.glob(f"{results_dir_protnormed}/*results.tsv")
    for file in all_files:
        df = pd.read_csv(file, sep="\t")
        assert ptm_normalization_tester.fcs_match_resultsconfig(list(df['log2fc']))


""" Check dampened ptms are as expected"""


@pytest.fixture
def data_paths():
    return {
        'changed': {
            'ptm': "../../test_data/unit_tests/protein_ptm_normalization/toy_data/results_ptm/condpair1.results.tsv",
            'protnormed': "../../test_data/unit_tests/protein_ptm_normalization/toy_data/results_ptm_protnormed/condpair1.results.tsv"
        },
        'unchanged': {
            'ptm': "../../test_data/unit_tests/protein_ptm_normalization/toy_data/results_ptm/condpair30.results.tsv",
            'protnormed': "../../test_data/unit_tests/protein_ptm_normalization/toy_data/results_ptm_protnormed/condpair30.results.tsv"
        }
    }

def compare_p_values(ptm_path, protnormed_path):
    results_df_ptm = pd.read_csv(ptm_path, sep="\t")
    results_df_protnormed = pd.read_csv(protnormed_path, sep="\t")
    
    ptm_p_values = np.log10(results_df_ptm["p_value"])
    protnormed_p_values = np.log10(results_df_protnormed["p_value"])
    
    return ptm_p_values, protnormed_p_values

@pytest.mark.parametrize("condition", ['changed', 'unchanged'])
def test_ptm_normalization(data_paths, condition):
    paths = data_paths[condition]
    ptm_p_values, protnormed_p_values = compare_p_values(paths['ptm'], paths['protnormed'])
    
    assert min(ptm_p_values) <= min(protnormed_p_values), f"Minimum p-value condition not met for {condition} condition"
    
    # You can add more assertions here if needed, for example:
    assert np.mean(ptm_p_values) <= np.mean(protnormed_p_values), f"Mean p-value condition not met for {condition} condition"