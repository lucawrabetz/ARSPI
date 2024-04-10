import os
import pdb
import time
import sys
import pandas as pd
import networkx as nx
from collections import defaultdict
from itertools import combinations

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
INSTANCES = '/home/luw28/Dropbox/ARSPI/exp/instances/'
RESULTS_DIR = 'results'

class CompoundFeatureComputer:
    """
    Object to  
    Args:
        row (pd series, row of the experiment df)
    Attributes:
        nodes, arcs, source, sink, policies, followers (int) - basic attributes
        row (series)
        set_name (string) - base set name
        instance_name (string) - base instance name
        graph (string) - path to graph edge list file
        cost (string) - path to arc costs file
        linked_list_arc_indices (dict) - only populated if self.construct_linkedlists_edges is called
        rev_linked_list_arc_indices (dict) - only populated if self.construct_linkedlists_edges is called
        clusters (list) - only populated if self.construct_clusters is called
    """
    def __init__(self, row):
        self.instance_set_path = os.path.join(INSTANCES, row["set_name"])
        self.nodes = row['nodes']
        self.arcs = row['arcs']
        self.source = 0
        self.sink = self.nodes - 1
        self.row = row
        self.policies = row['policies']
        self.followers = row['scenarios']
        self.set_name = row['set_name']
        self.instance_name = row['instance_name']
        self.graph = os.path.join(self.instance_set_path, self.instance_name.split(
            '-')[0] + '-' + self.instance_name.split('-')[1] + ".txt")
        self.costs = os.path.join(self.instance_set_path, self.instance_name.split('-')[0] + '-' + self.instance_name.split(
            '-')[1] + '-costs_' + self.instance_name.split('-')[2] + '.csv')
        self.linked_list_arc_indices = defaultdict(dict)
        self.rev_linked_list_arc_indices = defaultdict(dict)
        self.contracted_linked_list_arc_indices = defaultdict(dict)
        self.contracted_rev_linked_list_arc_indices = defaultdict(dict)
        self.clusters = [[] for _ in range(self.policies)]
        self.all_paths_total_costs = []


    def construct_clusters(self):
        # TODO (lucawrabetz): need to find a place to check or handle case where policies = 0 (currently checked at call site)
        assignments = self.row['partition'].split('-')
        follower = 0
        for a in assignments:
            self.clusters[int(a)].append(follower)
            follower += 1

    def construct_linkedlists_edges(self):
        self.G = nx.DiGraph()
        self.G_st = nx.DiGraph()
        self.contracted_node = self.source
        with open(self.graph, "r") as file:
            a = 0
            for line in file:
                line = line.strip()
                u, v = map(int, line.split())
                # assume no edges s-t exist, assume no edges such that u = t, or v = s exist
                # always add the edge to the normal graph
                self.G.add_edge(u, v)
                if u != self.source and v != self.sink:
                    # 'normal' edge, add u, v to contracted graph
                    self.G_st.add_edge(u, v)
                    self.contracted_linked_list_arc_indices[u][v] = a
                    self.contracted_rev_linked_list_arc_indices[v][u] = a
                if u == self.source:
                    # source, v node - add contracted_node, v to contracted graph
                    self.G_st.add_edge(self.contracted_node, v)
                    self.contracted_linked_list_arc_indices[self.contracted_node][v] = a
                    self.contracted_rev_linked_list_arc_indices[v][self.contracted_node] = a
                if v == self.sink:
                    # u, sink node - add u, contracted_node to contracted graph
                    self.G_st.add_edge(u, self.contracted_node)
                    self.contracted_linked_list_arc_indices[u][self.contracted_node] = a
                    self.contracted_rev_linked_list_arc_indices[self.contracted_node][u] = a
                self.linked_list_arc_indices[u][v] = a
                self.rev_linked_list_arc_indices[v][u] = a
                a += 1

    def compute_all_paths_total_costs(self):
        df = pd.read_csv(self.costs, header=None)
        followers = df.shape[0]
        G = nx.DiGraph()
        linked_list_arc_indices = defaultdict(dict)
        with open(self.graph, "r") as file:
            # Loop over each line in the file
            a = 0
            for line in file:
                # Strip the trailing newline character
                line = line.strip()
                # Split the line by space to get the values of u and v
                u, v = map(int, line.split())
                G.add_edge(u, v)
                linked_list_arc_indices[u][v] = a
                a += 1
        source = 0
        sink = max(G.nodes())
        all_paths_total_costs = []
        print('starting paths for graph: ', self.graph)
        for path in nx.all_simple_paths(G, source, sink):
            path_total_costs = []
            path_arc_indices = []
            for i in range(len(path)-1):
                u = path[i]
                v = path[i+1]
                path_arc_indices.append(linked_list_arc_indices[u][v])
            for q in range(followers):
                cost = 0
                for a in path_arc_indices:
                    cost += df.iloc[q, a]
                path_total_costs.append(cost)
            all_paths_total_costs.append(path_total_costs)
        print('done with paths for graph: ', self.graph)
        self.all_paths_total_costs = all_paths_total_costs

    def compute_exact_alpha(self):
        self.construct_clusters()
        self.compute_all_paths_total_costs()
        alpha = sys.maxsize
        for c in self.clusters:
            for i in range(len(c)):
                for j in range(i+1, len(c)):
                    for path in self.all_paths_total_costs:
                        val = min(path[i], path[j]) / max(path[i], path[j])
                        if val < alpha:
                            alpha = val
        return alpha

    def all_Gstpaths_oflength_kappa(self, kappa):
        pass

    def compute_alphahat_kappa(self, kappa=3):
        # OPTION FOR COMPUTING FOR A GENERAL KAPPA
        costs_df = pd.read_csv(self.costs, header=None)
        self.construct_clusters()
        self.construct_linkedlists_edges()
        alpha = sys.maxsize
        pdb.set_trace()
        # paths_st_kappa = all_Gstpaths_oflength_kappa(kappa)
        paths_st_kappa = [
            [0, 1],
            [1, 2],
            [2, 3]
        ]
        for path in paths_st_kappa:
            # paths_st_kappa is expected to contain paths expressed in terms of arc indices
            for c in self.clusters:
                if len(c) == 1:
                    continue
                for i, j in combinations(c, 2):
                    c_i = 0
                    c_j = 0
                    for a in path:
                        c_i += costs_df.iloc[i, a]
                        c_j += costs_df.iloc[j, a]
                    val = min(c_i, c_j) / max(c_i, c_j)
                    if val < alpha:
                        alpha = val
        return alpha

    def populate_datastructures(self):
        self.construct_clusters()
        self.construct_linkedlists_edges()

    def compute_alphahat1(self):
        alpha = 1
        costs_df = pd.read_csv(self.costs, header=None)
        self.construct_clusters()
        follower = 0
        for c in self.clusters:
            if len(c) == 1:
                continue
            for i, j in combinations(c, 2):
                ratios = costs_df.apply(lambda col: min(
                    col[i], col[j]) / max(col[i], col[j]))
                val = ratios.min()
                if val < alpha:
                    alpha = val
        return alpha

    def compute_alphahat2(self):
        alpha = 1
        costs_df = pd.read_csv(self.costs, header=None)
        self.construct_clusters()
        self.construct_linkedlists_edges()
        for c in self.clusters:
            if len(c) == 1:
                continue
            for i, j in combinations(c, 2):
                for v in range(self.nodes):
                    # for every node v, we enumerate the pairs of arcs (u,v), (v, ell)
                    for a1 in self.rev_linked_list_arc_indices[v].values():
                        for a2 in self.linked_list_arc_indices[v].values():
                            num = min(
                                costs_df.iloc[i, a1] + costs_df.iloc[i, a2], costs_df.iloc[j, a1] + costs_df.iloc[j, a2])
                            den = max(
                                costs_df.iloc[i, a1] + costs_df.iloc[i, a2], costs_df.iloc[j, a1] + costs_df.iloc[j, a2])
                            val = num / den
                            if val < alpha:
                                alpha = val
                # we enumerate the pairs of arcs (s, u), (ell, t), where s and t are the source and sink nodes
                for a1 in self.linked_list_arc_indices[self.source].values():
                    for a2 in self.rev_linked_list_arc_indices[self.sink].values():
                        num = min(costs_df.iloc[i, a1] + costs_df.iloc[i, a2],
                                  costs_df.iloc[j, a1] + costs_df.iloc[j, a2])
                        den = max(costs_df.iloc[i, a1] + costs_df.iloc[i, a2],
                                  costs_df.iloc[j, a1] + costs_df.iloc[j, a2])
                        val = num / den
                        if val < alpha:
                            alpha = val
        return alpha

    def compute_alphahat(self, kappa: int = 1):
        alpha = 1
        costs_df = pd.read_csv(self.costs, header=None)
        self.populate_datastructures()
        for c in self.clusters:
            if len(c) == 1:
                continue
            for i, j in combinations(c, 2):
                ratios = costs_df.apply(lambda col: min(
                    col[i], col[j]) / max(col[i], col[j]))
                val = ratios.min()
                if val < alpha:
                    alpha = val
        return alpha

def main():
    pass
    # results_file = BATCH_NAME + '.csv'
    # results_filepath = os.path.join(RESULTS_DIR, results_file)
    # results_df = pd.read_csv(results_filepath)
    # results_df = compute_approximation_ratio(
    #     results_df, exact="ENUMERATION_objective")
    # results_df = add_alphahat1(results_df)
    # results_df = add_alphahat2(results_df)
    # # results_df = add_alphahat_kappa(batch_exp_df, 2)
    # write_path = BATCH_NAME + '-alpha.csv'
    # results_df.to_csv(write_path)


if __name__ == '__main__':
    main()
