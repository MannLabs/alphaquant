from alphaquant.config.variables import QUANT_ID
import pandas as pd
import scipy.stats
import numpy as np
import sklearn.metrics

def find_node_parent_at_level(node, level):
    if node.level == level:
        return node
    while node.parent is not None:
        node = node.parent
        if node.level == level:
            return node

def find_node_parent_at_type(node, nodetype):
    if node.type == nodetype:
        return node
    while node.parent is not None:
        node = node.parent
        if node.type == nodetype:
            return node

def annotate_precursor_nodes(y_pred_total, ionnames_total, precursor_nodes):
    precursor2ml_score = {x:y for x, y in zip(ionnames_total, y_pred_total)}
    for precursor in precursor_nodes:
        ml_score = precursor2ml_score.get(precursor.name)
        precursor.ml_score = ml_score


def collect_node_parameters(all_nodes, w_annot = True):
    ion2param2val = {}
    all_headers = set()
    for node in all_nodes:

        param2val = get_param2val(node)
        if w_annot:
            param2val.update({QUANT_ID : node.name})
        if node.type == 'mod_seq_charge':
            if len(node.children) ==2:
                fcfc_diff = abs(node.children[0].fc - node.children[1].fc)
                param2val.update({"ms1_ms2_fcfc_diff" : fcfc_diff})

            for child in node.children:
                name_frac_mainclust = f"child_type_{child.type}_frac_mainclust"
                name_num_mainclusts = f"child_type_{child.type}_num_mainclusts"
                name_variance = f"child_type_{child.type}_cv_fcs"
                name_replicate_cv = f"child_type_{child.type}_replicate_cv"
                param2val.update({name_frac_mainclust: child.frac_mainclust, name_num_mainclusts : child.num_mainclusts, name_variance : calc_variance_for_node(child), name_replicate_cv : child.cv})
                clusterstats_dict = calc_cluster_stats(child)
                param2val.update(clusterstats_dict)
        all_headers.update(param2val.keys())
        ion2param2val[node.name] = param2val

    return get_dataframe(all_nodes,ion2param2val, all_headers)


def get_param2val(node):
    headers_of_interest = ["frac_mainclust", "num_mainclusts", "fraction_consistent",  "cv", "min_intensity", "min_reps"]
    param_dict = node.__dict__
    param2val = {x: param_dict.get(x) for x in headers_of_interest}
    param2val["num_leaves"] = len(node.leaves)
    if len(node.children)>0:
        param2val["num_children"] = len(node.children)
        param2val["cv_fcs"] = calc_variance_for_node(node)
    return param2val

def get_dataframe(all_nodes,ion2param2val, all_headers):
    all_headers = list(all_headers)
    rows = [[ion2param2val.get(node.name).get(header, np.nan) for header in all_headers] for node in all_nodes]
    df = pd.DataFrame(rows, columns = all_headers)
    return df


def calc_cluster_stats(node):
    num_elems_secondclust = 0
    fcs_clust0 = []
    fcs_clust1 = []
    intensities = []
    for child in node.children:
        if child.cluster ==0:
            fcs_clust0.append(child.fc)
        if child.cluster ==1:
            num_elems_secondclust+=1
            fcs_clust1.append(child.fc)
        intensities.append(child.min_intensity)

    if len(fcs_clust1)>0:
        betweenclust_fcfc = abs(np.mean(fcs_clust0) - np.mean(fcs_clust1))
    else:
        betweenclust_fcfc = 8

    stats_dict = {f"child_type_{node.type}_num_elems_secondclust" : num_elems_secondclust, f"child_type_{node.type}_num_clusters_total" : node.num_clusters, f"child_type_{node.type}_betweenclust_fcfc" : betweenclust_fcfc,
    f"child_type_{node.type}_num_elems_secondclust" : np.mean(intensities)}
    return stats_dict

def calc_variance_for_node(node):

    fcs_children = [x.fc for x in node.children]
    min_fc = min(fcs_children)
    fcs_children = [x.fc - min_fc for x in node.children]
    #LOGGER.info(f"fcs children are {fcs_children}, variance is {np.var(fcs_children)}")

    return scipy.stats.variation(fcs_children)
    

def calculate_average_r2_score(test_set_predictions):
    r2_scores = [sklearn.metrics.r2_score(y_true, y_pred) for y_true, y_pred in test_set_predictions]
    return np.mean(r2_scores)

