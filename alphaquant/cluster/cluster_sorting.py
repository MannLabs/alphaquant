
def decide_cluster_order(childnode2clust_init):
    """ranks the clusters from 0 to n (with 0 being the best) depending on the properties/similarities of the child nodes contained in each cluster
    """
    clust2score, clust2childnodes = get_clust2score_clust2childnodes(childnode2clust_init) #maps the idx of the cluster to the score and the childnodes of the cluster
    sorted_childnodes = get_sorted_childnodes_by_score(clust2score, clust2childnodes) #sorts the grouped childnodes by score
    childnode2clust = get_childnode2clust(sorted_childnodes) #maps the childnodes to the new cluster idx

    return childnode2clust

def get_clust2score_clust2childnodes(childnode2clust):
    clust2score = {}
    clust2childnodes = {}
    for childnode,clust in childnode2clust:
        clust2score[clust] = clust2score.get(clust, 0) + calculate_score(childnode)
        clust2childnodes[clust] = clust2childnodes.get(clust, []) + [childnode]
    return clust2score, clust2childnodes


def calculate_score(childnode):
    if childnode.type == "precursor": #boost the score of the cluster if it contains the precursor
        return 100
    return childnode.fraction_consistent*len(childnode.leaves)



def get_sorted_childnodes_by_score(clust2score, clust2childnodes, sort_descending = True): #higher score is better
    childnodes_score = [(childnodes, clust2score.get(clust)) for clust, childnodes in clust2childnodes.items()]
    childnodes_score = sorted(childnodes_score, key= lambda x : x[1], reverse=sort_descending)
    sorted_childnodes = [childnodes_score[i][0] for i in range(len(childnodes_score))]
    return sorted_childnodes

def get_childnode2clust(sorted_childnodes):
    childnode2clust = {}
    for clust_idx, childnodes in enumerate(sorted_childnodes):
        for childnode in childnodes:
            childnode2clust[childnode] = clust_idx
    return childnode2clust

