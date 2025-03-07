

class QuantifiedMultiConditionComparison():
    """
    class that contains all information to be used for the multi-condition
    analysis.
    """
    def __init__(self):
        self._proteinname2multicondprots = {}

    def add_quantified_proteins_to_comparison(self, quantified_proteins):
        for quantified_protein in quantified_proteins:
            self._initialize_multicond_protein_if_necessary(quantified_protein)
            self._add_quantified_protein_to_multicond(quantified_protein)

    def get_quantified_protein(self, protein_name, condpair):
        multicondprot =  self._proteinname2multicondprots.get(protein_name)
        return multicondprot.condpair2quantified_protein.get(condpair)


    def _initialize_multicond_protein_if_necessary(self, quantified_protein):
        if quantified_protein.name not in self._proteinname2multicondprots:
            self._proteinname2multicondprots[quantified_protein.name] = QuantifiedProteinMultiCondition(quantified_protein.name)

    def _add_quantified_protein_to_multicond(self, quantified_protein):
        multicondprot = self._proteinname2multicondprots.get(quantified_protein.name)
        multicondprot.condpair2quantified_protein[quantified_protein.condpair] = quantified_protein

class QuantifiedProteinMultiCondition():
    def __init__(self, name):
        self.name = name
        self.condpair2quantified_protein = {}


class QuantifiedProteinInCondpair():
    def __init__(self, condpair, name, log2fc, p_value, fdr):
        self.condpair = condpair
        self.name = name
        self.log2fc = log2fc
        self.p_value = p_value
        self.fdr = fdr


class QuantifiedProteinCanditateInCondpair(QuantifiedProteinInCondpair):
    """stump for cluster-based protein processing"""
    def __init__(self, peptide_nodes_of_cluster):
        self.peptides = None

    def _get_quantprot_properties_from_peptide_nodes(self, peptide_nodes_of_cluster):
        pass


# Cell

class ResultsDirectoryReader():
    def __init__(self, results_dir, condpairs_selected):
        self.quantified_multicondition_comparison = QuantifiedMultiConditionComparison() #initialize empty multicomparison object

        self.__condpairs_selected = condpairs_selected
        self.__localizer = ResultstableLocalizer(results_dir)
        self._add_all_condpairs_to_multicondition_comparison()

    def _add_all_condpairs_to_multicondition_comparison(self):
        for condpair in self.__condpairs_selected:
            file = self.__localizer.condpairname2file.get(condpair)
            quantified_proteins = ResultsTableReader(condpair, file).quantified_proteins
            self.quantified_multicondition_comparison.add_quantified_proteins_to_comparison(quantified_proteins)

# Cell
import pandas as pd

class ResultsTableReader():
    def __init__(self, condpair, file):
        self.quantified_proteins = []

        self.__condpair = condpair
        self.__results_df = self._read_table(file)
        self.__protnames = self._get_property("protein")
        self.__log2fcs = self._get_property("log2fc")
        self.__p_values = self._get_property("pval")
        self.__fdrs = self._get_property("fdr")

        self._init_quantified_proteins()

    def _read_table(self, file):
        return pd.read_csv(file, sep = "\t")

    def _get_property(self, property):
        return list(self.__results_df[property])

    def _init_quantified_proteins(self):
        for idx in range(len(self.__protnames)):
            self.quantified_proteins.append(QuantifiedProteinInCondpair(condpair=self.__condpair, name = self.__protnames[idx],
            log2fc= self.__log2fcs[idx], p_value=self.__p_values[idx], fdr = self.__fdrs[idx]))





# Cell
import re
class ResultstableLocalizer():
    def __init__(self, results_dir):
        self._resultsfilepaths = ResultsFiles(results_dir=results_dir).filepaths
        self.condpairname2file = {}
        self._load_condpairname2file()

    def _load_condpairname2file(self):
        for filepath in self._resultsfilepaths:
            filename = self._parse_condpairname_from_filepath(filepath)
            self.condpairname2file.update({filename: filepath})

    @staticmethod
    def _parse_condpairname_from_filepath(file):
        pattern = "(.*\/|^)(results.*\/)(.*)(.results.tsv)"
        matched = re.search(pattern, file)
        return matched.group(3)


import glob
class ResultsFiles():
    def __init__(self, results_dir):
        self._results_dir = results_dir
        self.filepaths = self.__get_result_filepaths()

    def __get_result_filepaths(self):
        return glob.glob(f'{self._results_dir}/*.results.tsv')

