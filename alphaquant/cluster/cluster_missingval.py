import alphaquant.cluster.cluster_ions as aq_cluster_ions
import anytree
import alphaquant.cluster.cluster_utils as aq_cluster_utils
import numpy as np
import scipy.stats
import numpy as np
import statistics


def create_protnode_from_missingval_ions(gene_name,diffions, normed_c1, normed_c2):
    return MissingValProtNodeCreator(gene_name, diffions, normed_c1, normed_c2).prot_node

class MissingValProtNodeCreator:
    def __init__(self, gene_name, diffions, normed_c1, normed_c2):
        self._gene_name = gene_name
        self._diffions = diffions
        self._normed_c1 = normed_c1
        self._normed_c2 = normed_c2

        self._nrep_c1 = None
        self._nrep_c2 = None
        self._all_intensities_c1 = None
        self._all_intensities_c2 = None
        self._total_intensity = None

        self._define_condition_properties()
        self.prot_node = self._create_protnode_from_missingval_ions()

    def _define_condition_properties(self):
        self._nrep_c1 = self._normed_c1.num_replicates
        self._nrep_c2 = self._normed_c2.num_replicates
        self._all_intensities_c1 = self._normed_c1.all_intensities
        self._all_intensities_c2 = self._normed_c2.all_intensities
        self._total_intensity = (np.mean(self._all_intensities_c1) +np.mean(self._all_intensities_c2))/2
    

    def _create_protnode_from_missingval_ions(self):
        #nrep_c1 and nrep_c2 are the number of replicates in the conditions in general, not the minimum required


        root_node = aq_cluster_ions.create_hierarchical_ion_grouping(self._gene_name, self._diffions)
        aq_cluster_ions.add_reduced_names_to_root(root_node)


        self._assign_properties_to_missingval_base_ions(root_node)
        nodes_to_test = self._get_nodes_to_test(root_node)
        levelname_nodes_to_test = nodes_to_test[0].level
        self._propagate_properties_to_nodes_to_test(nodes_to_test)
        self._propagate_properties_from_nodes_to_test_to_root(root_node, levelname_nodes_to_test)

        return root_node


    def _assign_properties_to_missingval_base_ions(self, root_node):

        for leaf in root_node.leaves:
            log2intensities_c1 = self._normed_c1.ion2nonNanvals.get(leaf.name)
            log2intensities_c2 = self._normed_c2.ion2nonNanvals.get(leaf.name)
            leaf.numvals_c1 = len(log2intensities_c1)
            leaf.numvals_c2 = len(log2intensities_c2)

            leaf.fc = np.nan
            
            leaf.missingval = True
            leaf.total_intensity = self._total_intensity
            leaf.fraction_consistent = np.nan

            leaf.fraction_missingval_c1 = self._normed_c1.ion2background.get(leaf.name).fraction_missingval
            leaf.fraction_missingval_c2 = self._normed_c2.ion2background.get(leaf.name).fraction_missingval

            leaf.intensity_c1 = np.mean(log2intensities_c1)
            leaf.intensity_c2 = np.mean(log2intensities_c2)


    @staticmethod
    def _get_nodes_to_test(root_node): #get the nodes in the lowest level that is relevant for the binomial test
        if root_node.depth == 2:
            return root_node.children
        else:
            return anytree.search.findall(root_node, filter_=lambda node: node.type == "mod_seq_charge")


    def _propagate_properties_to_nodes_to_test(self,nodes_to_test): #goes through each node to test and merges the properties from it's base to the node itself
        for node in nodes_to_test:
            for level_nodes in aq_cluster_utils.iterate_through_tree_levels_bottom_to_top(node):
                if level_nodes[0].level == "base":
                    continue

                for level_node in level_nodes:
                    self._aggregate_node_properties_missingval(level_node)
                    level_node.missingval = True




    def _propagate_properties_from_nodes_to_test_to_root(self, root_node, levelname_nodes_to_test):
        level_above_nodes_to_test = False
        for level_nodes in aq_cluster_utils.iterate_through_tree_levels_bottom_to_top(root_node):
            if level_nodes[0].level == levelname_nodes_to_test:
                self._assign_missingvals_prob_per_node(level_nodes)
                level_above_nodes_to_test = True
                continue
            if level_above_nodes_to_test:
                for level_node in level_nodes:
                    self._aggregate_node_properties_missingval(level_node)


    def _aggregate_node_properties_missingval(self, node):
        childs = node.children
        node.numvals_c1 = np.mean([child.numvals_c1 for child in childs])
        node.numvals_c2 = np.mean([child.numvals_c2 for child in childs])
        node.fraction_missingval_c1 = np.mean([child.fraction_missingval_c1 for child in childs])
        node.fraction_missingval_c2 = np.mean([child.fraction_missingval_c2 for child in childs])
        node.fc = np.mean([child.fc for child in childs])
        node.missingval = True
        node.fraction_consistent = np.nan
        node.total_intensity = np.sum([child.total_intensity for child in childs])
        node.intensity_c1 = np.mean([child.intensity_c1 for child in childs])
        node.intensity_c2 = np.mean([child.intensity_c2 for child in childs])
        if hasattr(childs[0], "z_val"):
            node.z_val = aq_cluster_utils.sum_and_re_scale_zvalues([child.z_val for child in childs])
            node.p_val = aq_cluster_utils.transform_znormed_to_pval(node.z_val)


    def _assign_missingvals_prob_per_node(self, nodes_to_test):
        for node in nodes_to_test:
            missingval_node_tester = MissingValNodeTester(node, self._nrep_c1, self._nrep_c2, self._all_intensities_c1, self._all_intensities_c2)
            node.p_val = missingval_node_tester.pval
            node.fc = missingval_node_tester.fc
            flipped_pval = 1-0.5*node.p_val #the flipped pval is always larger than 0.5 and the closer to 1 is gets, the closer it goes to 0.5, while the smaller it gets, the closer it goes to 1. When we express this with the standard normal distribution, we are always on the right side of the distribution, so we can use the inv_cdf function to get a positive z-value equivalent to the p-value
            node.z_val = abs(statistics.NormalDist().inv_cdf(flipped_pval))
            #the p-value can be obtained again by applying the transformation: statistics.NormalDist().cdf(z)*2 - 1



class MissingValNodeTester:
    def __init__(self, node_to_test, nrep_c1, nrep_c2, all_intensities_c1, all_intensities_c2):
        self._numvals_higher_condition = None
        self._numvals_lower_condition = None
        self._fraction_missingval_higher_condition = None
        self._fraction_missingval_lower_condition = None
        self._nrep_higher_condition = None
        self._nrep_lower_condition = None
        self._estimated_intensity_lower_condition = None
        self._all_intensities_higher_condition = None

        self.pval = None
        self.fc = None

        self._define_higher_and_lower_condition(node_to_test, nrep_c1, nrep_c2, all_intensities_c1, all_intensities_c2)
        self._define_pvalue_by_iterative_testing()
        self._define_matching_fc(node_to_test)

    
    def _define_higher_and_lower_condition(self, node_to_test, nrep_c1, nrep_c2, all_intensities_c1, all_intensities_c2):
        if node_to_test.numvals_c1 > node_to_test.numvals_c2:
            self._numvals_higher_condition = node_to_test.numvals_c1
            self._numvals_lower_condition = node_to_test.numvals_c2
            self._fraction_missingval_higher_condition = node_to_test.fraction_missingval_c1
            self._fraction_missingval_lower_condition = node_to_test.fraction_missingval_c2
            self._nrep_higher_condition = nrep_c1
            self._nrep_lower_condition = nrep_c2
            self._all_intensities_higher_condition = all_intensities_c1
        
        elif node_to_test.numvals_c2 > node_to_test.numvals_c1:
            self._numvals_higher_condition = node_to_test.numvals_c2
            self._numvals_lower_condition = node_to_test.numvals_c1
            self._fraction_missingval_higher_condition = node_to_test.fraction_missingval_c2
            self._fraction_missingval_lower_condition = node_to_test.fraction_missingval_c1
            self._nrep_higher_condition = nrep_c2
            self._nrep_lower_condition = nrep_c1
            self._all_intensities_higher_condition = all_intensities_c2

        else:
            raise Exception("Condition 1 and condition 2 have the same number of values. This should not be handled by the counting statistics module.")
    
    def _define_pvalue_by_iterative_testing(self):
        if self._perform_binomal_test_on_higher_condition() > 0.1: #the function returns a p-value
            self.pval = self._perform_binomal_test_on_lower_condition()
        
        else:
            self.pval = self._perform_fishers_exact_test()

    def _perform_binomal_test_on_higher_condition(self): # we first test the null hypothesis that the values observed in the higher condition (e.g. 5 values are there and we have 6 measurements in total) are missing at random. If this is not the case, we can't apply the binomial test to the lower condition. 
            pval_higher_condition = scipy.stats.binomtest(int(self._numvals_higher_condition), self._nrep_higher_condition, 1-self._fraction_missingval_higher_condition).pvalue
            return pval_higher_condition
    
    def _perform_binomal_test_on_lower_condition(self):
        pval_lower_condition = scipy.stats.binomtest(int(self._numvals_lower_condition), self._nrep_lower_condition, 1-self._fraction_missingval_higher_condition).pvalue
        return pval_lower_condition


    def _perform_fishers_exact_test(self):
        num_missing_higher_condition = self._nrep_higher_condition - self._numvals_higher_condition
        num_missing_lower_condition = self._nrep_lower_condition - self._numvals_lower_condition

        contingency_table = np.array([[self._numvals_higher_condition, num_missing_higher_condition],
                                    [self._numvals_lower_condition, num_missing_lower_condition]])
        
        odds_ratio, p = scipy.stats.fisher_exact(contingency_table)

        return p


    def _define_matching_fc(self, node_to_test):
        intensity_lower = np.percentile(self._all_intensities_higher_condition, 0.1) #simply sample from the lower 10% of intensities from the higher condition
        numvals_c1 = node_to_test.numvals_c1
        numvals_c2 = node_to_test.numvals_c2
        if numvals_c1 > numvals_c2:
            self.fc = node_to_test.intensity_c1 - intensity_lower
        elif numvals_c2 > numvals_c1:
            self.fc = intensity_lower - node_to_test.intensity_c2
        else:
            raise Exception("Condition 1 and condition 2 have the same number of values. This should not be handled by the binomial test.")
        

        