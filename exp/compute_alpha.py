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

# def alphahat_for_two_followers(costs_df, f1, f2):
#     ratios = costs_df.apply(lambda col: min(
#         col[f1], col[f2]) / max(col[f1], col[f2]))
#     # if ratios.min() < 0.5:
#     #     print(ratios)
#     return ratios.min()
# 
# 
# def compute_alphahat(row, instance_directory):
#     clusters = [[] for i in range(row['policies'])]
#     assignments = row['GREEDY_partition'].split('-')
#     follower = 0
#     instance_name = row['instance_name']
#     instance_file = instance_name.split('-')[0] + '-' + instance_name.split(
#         '-')[1] + '-costs_' + instance_name.split('-')[2] + '.csv'
#     costs_df = pd.read_csv(os.path.join(
#         instance_directory, instance_file), header=None)
#     alpha = sys.maxsize
#     for a in assignments:
#         clusters[int(a)].append(follower)
#         follower += 1
#     for c in clusters:
#         for i in range(len(c)):
#             for j in range(i+1, len(c)):
#                 val = alphahat_for_two_followers(costs_df, c[i], c[j])
#                 if val < alpha:
#                     alpha = val
#     return alpha
# 
# 
# def add_alphahat(df, exp_name):
#     print("adding alphahat 1...")
#     instance_directory = os.path.join('instances', exp_name)
#     # Iterate through each row and compute the value for the new column
#     new_column = []
#     times_column = []
#     for index, row in df.iterrows():
#         begin_seconds = time.time()
#         new_value = compute_alphahat(row, instance_directory)
#         duration_seconds = time.time() - begin_seconds
#         new_column.append(new_value)
#         times_column.append(duration_seconds)
#         print("computed alphahat1 value ", new_value,
#               " in ", duration_seconds, " seconds.")
#     # Add the new column to the DataFrame
#     df['alpha_hat_1'] = new_column
#     df['alpha_hat_1_time'] = times_column
#     return df
# 
# 
# def construct_linkedlists_edges(instance_directory, instance_file, graph_file, P):
#     G = nx.DiGraph()
#     linked_list_arc_indices = defaultdict(dict)
#     rev_linked_list_arc_indices = defaultdict(dict)
#     with open(os.path.join(instance_directory, graph_file), "r") as file:
#         # Loop over each line in the file
#         a = 0
#         for line in file:
#             # Strip the trailing newline character
#             line = line.strip()
#             # Split the line by space to get the values of u and v
#             u, v = map(int, line.split())
#             G.add_edge(u, v)
#             linked_list_arc_indices[u][v] = a
#             rev_linked_list_arc_indices[v][u] = a
#             a += 1
#     return linked_list_arc_indices, rev_linked_list_arc_indices, max(G.nodes())+1, 0, max(G.nodes())
# 
# 
# def file_names(row, instance_directory):
#     instance_name = row['instance_name']
#     instance_file = instance_name.split('-')[0] + '-' + instance_name.split(
#         '-')[1] + '-costs_' + instance_name.split('-')[2] + '.csv'
#     graph_file = instance_name.split(
#         '-')[0] + '-' + instance_name.split('-')[1] + ".txt"
#     return instance_name, instance_file, graph_file
# 

class InstanceOutputRow:
    """
    InstanceOutputRow object to hold instance name, cost file, graph file, and row
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
        self.clusters = [[] for i in range(self.policies)]
        self.all_paths_total_costs = []


    def construct_clusters(self):
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

    def compute_alphahat1(self):
        costs_df = pd.read_csv(self.costs, header=None)
        self.construct_clusters()
        follower = 0
        alpha = sys.maxsize
        for c in self.clusters:
            if len(c) == 1:
                continue
            for i in range(len(c)):
                for j in range(i+1, len(c)):
                    ratios = costs_df.apply(lambda col: min(
                        col[i], col[j]) / max(col[i], col[j]))
                    val = ratios.min()
                    if val < alpha:
                        alpha = val
        return alpha

    def compute_alphahat2(self):
        costs_df = pd.read_csv(self.costs, header=None)
        self.construct_clusters()
        self.construct_linkedlists_edges()
        alpha = sys.maxsize
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

    def compute_alphahat_kappaAAA(self, kappa):
        # OPTION FOR COMPUTING FOR A GENERAL KAPPA
        costs_df = pd.read_csv(self.costs, header=None)
        self.construct_clusters()
        self.construct_linkedlists_edges()
        alpha = sys.maxsize
        pdb.set_trace()
        for c in self.clusters:
            if len(c) == 1:
                continue
            for i, j in combinations(c, 2):
                # for path in paths_kappa(G_stcontracted)
                # Step 2: Find all simple paths of size N in the contracted graph
                all_paths = list(nx.all_simple_paths(
                    self.G_st, self.contracted_node, self.contracted_node))
                paths_of_size_kappa = [
                    path for path in all_paths if len(path) == kappa]
                for path in paths_of_size_kappa:
                    path_cost_i = 0
                    path_cost_j = 0
                    for index in range(len(path) - 1):
                        u = path[index]
                        v = path[index+1]
                        a = self.linked_list_arc_indices[u][v]
                        path_cost_i += costs_df.iloc[i, a]
                        path_cost_j += costs_df.iloc[j, a]
                        val = min(path_cost_i, path_cost_j) / \
                            max(path_cost_i, path_cost_j)
                        if val < alpha:
                            alpha = val
        return alpha

def add_alphahat1(df):
    print("adding alphahat 1...")
    # Iterate through each row and compute the value for the new column
    new_column = []
    times_column = []
    for index, row in df.iterrows():
        instance = InstanceOutputRow(row)
        begin_seconds = time.time()
        new_value = instance.compute_alphahat1()
        duration_seconds = time.time() - begin_seconds
        new_column.append(new_value)
        times_column.append(duration_seconds)
        print("computed alphahat1 value ", new_value,
              " in ", duration_seconds, " seconds.")
    # Add the new column to the DataFrame
    df['alpha_hat_1'] = new_column
    df['alpha_hat_1_time'] = times_column
    return df

def add_alphahat2(df):
    print("adding alphahat 2...")
    alphahatn_column = []
    times_column = []
    for index, row in df.iterrows():
        instance = InstanceOutputRow(row)
        begin_seconds = time.time()
        new_value = instance.compute_alphahat2()
        duration_seconds = time.time() - begin_seconds
        alphahatn_column.append(new_value)
        times_column.append(duration_seconds)
        print("computed alphahat2 value ", new_value,
              " in ", duration_seconds, " seconds.")
    df['alpha_hat_2'] = alphahatn_column
    df['alpha_hat_2_time'] = times_column
    return df


def add_alphahat_kappa(batch_exp_df, kappa=3):
    print("adding alphahat ", kappa, "...")
    alphahat_column = []
    times_column = []
    for index, row in batch_exp_df.iterrows():
        instance = InstanceOutputRow(row)
        begin_seconds = time.time()
        new_value = instance.compute_alphahat_kappa(kappa)
        duration_seconds = time.time() - begin_seconds
        alphahat_column.append(new_value)
        times_column.append(duration_seconds)
        print("computed alphahat", kappa, " value ",
              new_value, " in ", duration_seconds, " seconds.")
    column_name = "alpha_hat_" + str(kappa)
    time_column_name = column_name + "_time"
    batch_exp_df[column_name] = alphahat_column
    batch_exp_df[time_column_name] = times_column
    return batch_exp_df


def compute_all_paths_total_costs(row, instance_directory, graph_file, instance_file):
    df = pd.read_csv(os.path.join(
        instance_directory, instance_file), header=None)
    followers = df.shape[0]
    G = nx.DiGraph()
    linked_list_arc_indices = defaultdict(dict)
    with open(os.path.join(instance_directory, graph_file), "r") as file:
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
    print('starting paths for graph: ', graph_file)
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
    print('done with paths for graph: ', graph_file)

    return all_paths_total_costs


def compute_exact_alpha(row, all_paths_total_costs):
    clusters = [[] for i in range(row['policies'])]
    assignments = row['GREEDY_partition'].split('-')
    follower = 0
    alpha = sys.maxsize
    for a in assignments:
        clusters[int(a)].append(follower)
        follower += 1
    for c in clusters:
        for i in range(len(c)):
            for j in range(i+1, len(c)):
                for path in all_paths_total_costs:
                    val = min(path[i], path[j]) / max(path[i], path[j])
                    if val < alpha:
                        alpha = val
    return alpha


def add_exact_alpha(df, exp_name, batch_directory):
    new_column = []
    all_paths_total_costs_dict = {}
    for index, row in df.iterrows():
        instance_name = row['instance_name']
        graph_file = instance_name.split(
            '-')[0] + '-' + instance_name.split('-')[1] + ".txt"
        instance_file = instance_name.split('-')[0] + '-' + instance_name.split(
            '-')[1] + '-costs_' + instance_name.split('-')[2] + '.csv'
        path_costs = all_paths_total_costs_dict.setdefault(
            graph_file, compute_all_paths_total_costs(row, batch_directory, graph_file, instance_file))
        new_value = compute_exact_alpha(row, path_costs)
        new_column.append(new_value)

    # Add the new column to the DataFrame
    df['alpha'] = new_column
    return df


def compute_approximation_ratio(df, exact="MIP_objective"):
    df['approximation_ratio'] = df['GREEDY_objective'] / df[exact]
    print('done with approximation ratio')
    return df


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
