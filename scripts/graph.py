import os
import random
import math
import pandas as pd
import numpy as np
import networkx as nx
from datetime import date

def append_date(basename):
    """
    Append today's date to experiment name
    """
    today = date.today()
    date_str = today.strftime("%m_%d_%y")
    name = basename + "-" + date_str
    return name


def check_make_dir(path, i):
    """
    Recursively check if an experiment directory exists, or create one with the highest number
        - example - if "path" string is "/dat/experiments/test-01_29_22", and there already exist:
            - "/dat/experiments/test-01_29_22-0"
            - "/dat/experiments/test-01_29_22-1"
            - "/dat/experiments/test-01_29_22-2"
        we have to create the dir "/dat/experiments/test-01_29_22-3"
    """
    isdir_full = os.path.isdir(path + "-" + str(i))
    # if the directory exists, call on the next i
    if isdir_full:
        return check_make_dir(path, i + 1)
    # base case - create directory for given i (and return final path)
    else:
        os.mkdir(path + "-" + str(i))
        return path + "-" + str(i)


# def arcs_for_current_layer(current_layer, next_layer, p):
#     '''
#     - add arcs for this layer - outgoing
#     - random arcs to nodes in all 'future layers'
#     - 'next' is just all the remaining nodes in future layers
#     - p is the probability that the node i is connected to any j in future layers - to simplify we use it as a proportion and convert to arcs_per_node
#     '''
# 
#     new_arcs = []
#     arcs_per_node = math.floor(len(next_layer) * p)
# 
#     if arcs_per_node < 1:
#         arcs_per_node = 1
# 
#     for i in current_layer:
#         try:
#             sampled_nodes = random.sample(next_layer, arcs_per_node)
#         except ValueError:
#             print("SAMPLE EXCEPTION OCCURRED")
#             print("len next_layer: " + str(len(next_layer)))
#             print("arcs_per_node: " + str(arcs_per_node))
#             print("p: " + str(p))
#         for j in sampled_nodes:
#             new_arc = (i, j)
#             new_arcs.append(new_arc)
# 
#     return new_arcs


class CompleteLayerInput:
    def __init__(self,
                 num_layers=5,
                 num_per_layer=10,
                 followers=5,
                 follower_groups=3,
                 #budget=3,
                 high_mean=200,
                 low_mean=50,
                 standard_deviation=10,
                 set_name="graphs",
                 num_cost_files=30):
        '''
            - Assertions (input sanity checks):
                - num_layers >= budget;
                - num_per_layer >= follower_groups;
                - create a directory for the set_directory if it doesn't exist.
        '''
        # Make assertions / sanity checks.
        #assert num_layers >= budget
        self.set_directory = os.path.join("dat", set_name)
        if not os.path.isdir(self.set_directory):
            os.mkdir(self.set_directory)
        self.num_layers = num_layers
        self.num_per_layer = num_per_layer
        self.followers = followers
        self.follower_groups = follower_groups
        #self.budget = budget
        self.high_mean = high_mean
        self.low_mean = low_mean
        self.standard_deviation = standard_deviation
        self.set_name = set_name
        self.num_cost_files = num_cost_files


class CompleteLayerGraph:
    def __init__(self,
                 instance_input):
        '''
        assign basic parameters
        '''
        self.source = 0
        self.arc_list = []
        self.arc_costs = []
        self.num_layers = instance_input.num_layers
        self.num_per_layer = instance_input.num_per_layer
        self.nodes = 2 + self.num_layers * self.num_per_layer
        self.sink = self.nodes - 1
        self.followers = instance_input.followers
        self.follower_groups = instance_input.follower_groups
        # self.budget = instance_input.budget
        self.high_mean = instance_input.high_mean
        self.low_mean = instance_input.low_mean
        self.standard_deviation = instance_input.standard_deviation
        self.set_name = instance_input.set_name
        self.set_directory = instance_input.set_directory
        self.G = None
        self.full_graph_name = self.set_name + "-" + str(self.nodes) + "_" + str(self.follower_groups)
        self.num_cost_files = instance_input.num_cost_files
        self.arc_costs = [[[] for i in range(self.followers)] for j in range(self.num_cost_files)]

    def populate_arc_list(self):
        '''
        populate the arc list (fully connect the graph) - each node in every layer has an
        arc to every node in the next layer
        '''
        i = self.source
        j = 1
        while j <= self.num_per_layer:
            self.arc_list.append((i, j))
            j += 1
        for k in range(self.num_layers - 1):
            current_layer = list(range(i + 1, i + self.num_per_layer + 1))
            next_layer = list(range(i + self.num_per_layer + 1, i + 2*self.num_per_layer + 1))
            for i in current_layer:
                for j in next_layer:
                    self.arc_list.append((i, j))
        for i in range(self.sink - self.num_per_layer, self.sink):
            self.arc_list.append((i, self.sink))
        self.arcs = len(self.arc_list)
        self.G = nx.DiGraph(self.arc_list)


    def generate_costs(self):
        '''
        generate costs based on the following rule:
            - for every group of followers select a random choice of r_0 arcs (no overlap between groups)
            - these arcs will have costs drawn from N(low_mean, standard_deviation) for that group
            - the rest of the arcs will have costs drawn from N(high_mean, standard_deviation) for the group
        '''
        for cost_id in range(self.num_cost_files):
            cheap_matrix = [[0 for a in range(self.arcs)] for i in range(self.follower_groups)]
            # arcs_per_layer = self.num_per_layer ** 2
            # for layer in range(1, 1+self.budget):
            #     # we are considering the first (node 0) and middle layers
            #     # num_layers refers to the middle layers
            #     if layer == 0:
            #         first = 0
            #         last = self.num_per_layer
            #     elif layer == self.num_layers:
            #         first = self.arcs - self.num_per_layer
            #         last = self.arcs
            #     else:
            #         first = (layer - 1) * arcs_per_layer + self.num_per_layer
            #         last = first + arcs_per_layer
            #     arc_set = set(range(first, last))

            arc_set = set(range(self.source, self.arcs))
            for group in range(self.follower_groups):
                cheap_arcs = set(random.sample(arc_set, self.follower_groups))
                for arc in cheap_arcs:
                    cheap_matrix[group][arc] = 1
            group_index = 0
            for i in range(self.followers):
                for a in range(self.arcs):
                    if cheap_matrix[group_index][a]:
                        cost = int(np.round(np.random.normal(self.low_mean, self.standard_deviation)))
                    else:
                        cost = int(np.round(np.random.normal(self.high_mean, self.standard_deviation)))
                    self.arc_costs[cost_id][i].append(cost)
                if group_index == self.follower_groups - 1: group_index = 0
                else: group_index += 1

    def write_costs(self):
        '''
        write the generated costs to file
        '''
        for cost_id in range(self.num_cost_files):
            costs_file_name = self.full_graph_name + "-costs_" + str(self.followers) + "_" + str(cost_id) + ".csv"
            costs_file_path = os.path.join(self.set_directory, costs_file_name)
            costs_df = pd.DataFrame(self.arc_costs[cost_id])
            costs_df.to_csv(costs_file_path, header=False, index=False)

    def write_graph(self):
        '''
        write the edge list to a file
        '''
        file_name = self.full_graph_name + ".txt"
        file_path = os.path.join(self.set_directory, file_name)
        nx.write_edgelist(self.G, file_path, data=False)

    def populate_graph(self, write=True):
        '''
        do all work outside of constructor: populate arc_list and generate_costs
        if param write is true (true by default) write graph and costs to files
        '''
        self.populate_arc_list()
        self.generate_costs()
        if write:
            self.write_costs()
            self.write_graph()


class LayerGraph:
    s = 0
    num_layers = 0
    num_per_layer = 0
    p = 0
    arcs = []
    n = 0
    m = 0
    t = 0

    def __init__(self, num_layerss, num_per_layerr, pp):
        '''
        only generates graph topology (no cost information in this class)
        self.p is the probability that a node is connected to any given node ahead of it
        - we'll use it as a proportion to calculate arcs per node and sample a random set of that number
        '''

        self.num_layers = num_layerss
        self.num_per_layer = num_per_layerr
        self.p = pp
        self.n = 2 + self.num_layers*self.num_per_layer
        self.t = self.n-1

        current_layer = [self.s]
        next_layer = [i for i in range(self.n) if i not in current_layer]
        # import pdb
        # pdb.set_trace()

        for i in range(self.num_layers+1):
            new_arcs = arcs_for_current_layer(
                current_layer, next_layer, self.p)

            self.arcs.extend(new_arcs)
            node1_next_layer = current_layer[-1] + 1
            current_layer = []

            for i in range(self.num_per_layer):
                current_layer.append(node1_next_layer + i)

            node1_rest = current_layer[-1]
            next_layer = [i for i in range(self.n) if i > node1_rest]

        self.m = len(self.arcs)

    def printGraph(self, edge_list=True):
        '''
        - print out the graph
        '''
        print("n: " + str(self.n) + " - " + str([i for i in range(self.n)]))
        print("m: " + str(self.m))
        if edge_list:
            print("arcs: ")
            # import pdb; pdb.set_trace()

            for i in range(self.m):
                print("     " + str(self.arcs[i]))

    def checksNX(self, filename):
        '''
        - check for uniqueness of edges to avoid multigraph issue
        - convert graph to a networkx object
        - check anything you want - connectivity, parallel edges, etc
        - also writes graph to a file in standard edge list style for ya
        '''

        edge_set = set(self.arcs)
        self.arcs = list(edge_set)
        self.m = len(self.arcs)

        G = nx.DiGraph(self.arcs)
        path = nx.has_path(G, self.s, self.t)
        nx.write_edgelist(G, filename, data=False)

        return path


class ErdosRenyi:
    s = -1
    t = -1
    n = -1
    pr = -1
    density = -1
    m = -1
    G = -1

    def __init__(self, nodes, probability):
        '''
        Generate graph topology using networkx erdos-renyi graph
        '''
        self.n = nodes
        self.pr = probability
        self.s = 0
        self.t = self.n-1

        diameter = -1
        st_shortest_path = -1

        while st_shortest_path < (diameter) or st_shortest_path < 0 or diameter < 0:
            self.G = nx.erdos_renyi_graph(self.n, self.pr, directed=True)

            try:
                diameter = nx.diameter(self.G)
            except:
                if (not nx.is_strongly_connected(self.G)):
                    diameter = -1

            try:
                st_shortest_path = len(nx.shortest_path(self.G, source=self.s, target=self.t)) - 1
            except:
                st_shortest_path = -1

        edge_list = [e for e in self.G.edges]
        self.m = len(edge_list)
        self.density = nx.density(self.G)

        print("Erdos Renyi Graph Generated")

    def printGraph(self, edge_list=True):
        '''
        - print out the graph
        '''
        print("n: " + str(self.n) + " - " + str([i for i in range(self.n)]))
        print("m: " + str(self.m))
        if edge_list:
            print("arcs: ")
            # import pdb; pdb.set_trace()

            for i in range(self.m):
                print("     " + str(self.arcs[i]))

    def addEdgeAttr(self, sub):
        '''
        - add value of sub (subgraph index) to all edges
        - loop through edges and add sub
        '''
        for edge in self.G.edges():
            i = edge[0]
            j = edge[1]
            self.G.edges[i, j]["sub"] = sub


    def checksNX(self, filename):
        '''
        - check anything you want - connectivity, parallel edges, etc
        - also writes graph to a file in standard edge list style
        '''
        nx.write_edgelist(self.G, filename, data=False)


class CompoundGraph:
    s = -1
    t = -1
    n = -1
    density = -1
    m = -1
    G = -1

    def __init__(self, graphs):
        '''
        Join graphs at source and sink nodes
            - graphs - list of graph objects
        '''
        # assume all graphs in graphs have same n
        n_0 = graphs[0].n
        self.n = n_0 * len(graphs) + 2
        self.s = 0
        self.t = self.n - 1

        # create source and sink networkx graphs
        sourceG = nx.DiGraph()
        sourceG.add_node(0)
        sinkG = nx.DiGraph()
        sinkG.add_node(self.t)
        temp_graphs = [0 for G in graphs]

        # loop through every graph in graphs and perform disjoint unions sequentially
        for i in range(len(graphs)):
            graphs[i].addEdgeAttr(i+1)
            if i == 0:
                temp_graphs[0] = nx.disjoint_union(sourceG, graphs[0].G)
            else:
                temp_graphs[i] = nx.disjoint_union(temp_graphs[i-1], graphs[i].G)

        self.G = nx.disjoint_union(temp_graphs[-1], sinkG)

        # add the edges to connect and source and sink
        # if we want to add more edges between subgraphs we can do it here as well
        t_0 = graphs[0].t
        for i in range(len(graphs)):
            new_edges = [(self.s, n_0 * (i) + 1, {"sub": 0}), (t_0 + n_0 * (i) + 1, self.t, {"sub": 0})]
            self.G.add_edges_from(new_edges)

        # final class attributes
        edge_list = [e for e in self.G.edges]
        self.m = len(edge_list)
        self.density = nx.density(self.G)

    def writeGraph(self, filename):
        '''
        - check anything you want - connectivity, parallel edges, etc
        - also writes graph to a file in standard edge list style
        '''
        nx.write_edgelist(self.G, filename, data=True)


class PseudoTreeRawInput:
    def __init__(self,
                 num_layers=1,
                 num_per_layer=2,
                 followers=1,
                 follower_groups=1,
                 # budget=3,
                 high_mean=200,
                 low_mean=50,
                 standard_deviation=10,
                 interdiction_delta=50,
                 set_name="graphs",
                 set_directory="dat/graphs"):
        '''
            - Assertions (input sanity checks):
                - num_layers >= budget;
                - num_per_layer >= follower_groups;
                - create a directory for the set_directory if it doesn't exist.
        '''
        # Make assertions / sanity checks.
        # assert num_layers >= budget
        # assert num_per_layer >= follower_groups
        if not os.path.isdir(set_directory):
            os.mkdir(set_directory)
        self.num_layers = num_layers
        self.num_per_layer = num_per_layer
        self.followers = followers
        self.follower_groups = follower_groups
        # self.budget = budget
        self.high_mean = high_mean
        self.low_mean = low_mean
        self.standard_deviation = standard_deviation
        self.interdiction_delta = interdiction_delta
        self.set_name = set_name
        self.set_directory = set_directory

class PseudoTreeInput:
    def __init__(self,
                 followers=1,
                 follower_groups=3,
                 layers=3, # not including sink node, just the rary tree (including root)
                 # budget=3,
                 high_mean=200,
                 low_mean=50,
                 standard_deviation=10,
                 interdiction_delta=50,
                 set_name="graphs",
                 set_directory="dat/graphs"):
        '''
            - Assertions (input sanity checks):
                - create a directory for the set_directory if it doesn't exist.
        '''
        # Make assertions / sanity checks.
        # assert num_layers >= budget
        # assert num_per_layer >= follower_groups
        if not os.path.isdir(set_directory):
            os.mkdir(set_directory)
        self.followers = followers
        self.follower_groups = follower_groups
        self.layers = layers
        # self.budget = budget
        self.high_mean = high_mean
        self.low_mean = low_mean
        self.standard_deviation = standard_deviation
        self.interdiction_delta = interdiction_delta
        self.set_name = set_name
        self.set_directory = set_directory

class PseudoTreeRaw:
    '''
        This graph class is designed to meet the condition that G[N \ t] is a tree.
        To generate this graph, we'll build it path by path, adding paths to the source node that are always disjoint.
        In order to make the tree interesting, we will also let the paths split into new paths, without rejoining which would induce a cycle.
        In order to make the instance interesting, we will have as many paths / path groups starting from the source node, as follower groups.
        Each follower will have its costs drawn from the cheap distribution on one of the path groups, and expensive on the rest.
        SOME IMPORTANT THINGS TO NOTE:
            - Layers means something different, than in the layer graph, i.e. it means the level at which a split occurs.
                - e.g. an instance with 3 layers, has 3 "splits"
                - additionally, num_per_layer will be the number of nodes between every split, between s and the first split,
                and between the last split and t.
            - That said, the input structure is very consistent with CompleteLayerGraph, however the budget does not matter for this graph.
    '''
    def __init__(self,
                 instance_input):
        self.source = 0
        self.arc_list = []
        self.arc_costs = []
        self.num_layers = instance_input.num_layers
        self.num_per_layer = instance_input.num_per_layer
        self.nodes = 0 # worked out when nodes are added manually
        self.sink = 0 # also worked out when self.nodes is set
        self.followers = instance_input.followers
        self.follower_groups = instance_input.follower_groups
        # self.budget = instance_input.budget
        self.high_mean = instance_input.high_mean
        self.low_mean = instance_input.low_mean
        self.standard_deviation = instance_input.standard_deviation
        self.interdiction_delta = instance_input.interdiction_delta
        self.set_name = instance_input.set_name
        self.set_directory = instance_input.set_directory
        self.G = None
        self.full_graph_name = None # worked out when self.nodes is set
        self.arc_costs = [[] for i in range(self.followers+1)] # Interdiction_costs are in last list / row.

    def create_path_group(self, layer, root, split_n):
        '''
            - create a path group starting from node sink_node, which could be the sink node of the graph or a splitting node
            - create the paths recursively at splits, splitting into split_n paths
            - add arcs to arc list
            - if layer == num_layers+1, stop
        '''
        # append all arcs to arc list for root to root + self.num_per_layer (straight line path)
        for i in range(root, root+self.num_per_layer):
            j = i+1 # needs to be adjusted for the splits
            self.arc_list.append((i, j))
        if layer <= self.num_layers:
            create_path_group(self, layer+1, j, 1)
        else: return

    def populate_arc_list(self):
        '''
            - populate the arc list
            - we will do this by adding to the nodes (and creating arcs along the way) along the path groups (splitting at every layer)
            - we follow the flow: create a path group
            - initialize graph, as well as self.nodes, self.full_graph_name, and self.sink, once the work is done
        '''
        create_path_group(self, layer, self.source, 1)
        # self.nodes = 
        # self.full_graph_name = self.set_name + "-" + str(self.nodes) + "_" + str(self.follower_groups)
        # self.sink = self.nodes - 1


class PseudoTree:
    '''
        This graph class is designed to meet the condition that G[N \ t] is a tree.
        We construct our graph using the networkx built-in full_rary_tree function.
    '''
    def __init__(self,
                 instance_input):
        self.source = 0
        self.arc_list = []
        self.arc_costs = []
        self.num_layers = instance_input.layers # layers includes root node but not sink
        # every layer l (where 0 is the root node) has follower_groups ^ l nodes,
        # so n is the sum of a geometric series sum_{l=0}{layers-1} follower_groups^l (+ 1 because we will add the sink node)
        self.nodes = int(1 + (1-instance_input.follower_groups**(instance_input.layers)) / (1 - instance_input.follower_groups))
        self.sink = self.nodes - 1
        self.followers = instance_input.followers
        self.follower_groups = instance_input.follower_groups
        # self.budget = instance_input.budget
        self.high_mean = instance_input.high_mean
        self.low_mean = instance_input.low_mean
        self.standard_deviation = instance_input.standard_deviation
        self.interdiction_delta = instance_input.interdiction_delta
        self.set_name = instance_input.set_name
        self.set_directory = instance_input.set_directory
        self.G = nx.full_rary_tree(self.follower_groups, self.nodes-1)
        self.G.add_node(self.sink)
        # layer structure is as follows, where l denotes the layer, and l=0 is the layer of the root/source node:
        #   - l=0 includes node 0
        #   - l=1 includes nodes 1, ...., follower_groups
        #   - l=2 includes nodes follower_groups+1, ..., follower_groups+follower_groups^2
        #   - l=l includes nodes sum_{i=0}^{l-1} follower_groups^i, ..., sum_{i=1}^{l} follower_groups^i
        # to add edges from last layer to sink, we iterate over leaf nodes:
        l = self.num_layers-1 # our formula is 0-indexed
        start_node = int((1-self.follower_groups**(l)) / (1 - self.follower_groups)) # first node of last layer
        end_node = int((self.follower_groups - self.follower_groups**(l+1)) / (1 - self.follower_groups)) # last node of last layer
        for i in range(start_node, end_node+1):
            self.G.add_edge(i, self.sink)
        self.full_graph_name = self.set_name + "-" + str(self.nodes) + "_" + str(self.follower_groups)
        self.arc_costs = [[] for i in range(self.followers+1)] # Interdiction_costs are in last list / row.


def compound_generator(n_0, pr, subgraphs):
    '''
    - helper function to generate a compound graph and return it
    - just generating the ErdosRenyi graphs for CompoundGraph constructor
    - n_0 and pr are the node number and probability value for the individual subgraphs
    '''
    graphs = [ErdosRenyi(n_0, pr) for i in range(subgraphs)]

    G = CompoundGraph(graphs)

    return G


def complete_layer_inputs():
    inputs = []
    inputs.append(CompleteLayerInput())
    inputs.append(CompleteLayerInput(
                 num_layers=10,
                 num_per_layer=4,
                 followers=10,
                 follower_groups=10,
                 #budget=5,
                 high_mean=200,
                 low_mean=50,
                 standard_deviation=10,
                 set_name="small_testbed_NE"))
    inputs.append(CompleteLayerInput(
                 num_layers=8,
                 num_per_layer=4,
                 followers=10,
                 follower_groups=10,
                 #budget=5,
                 high_mean=200,
                 low_mean=50,
                 standard_deviation=10,
                 set_name="tiny_testbed_NE"))
    return inputs

if __name__ == "__main__":
    for i in range(2, 3):
        inputs = complete_layer_inputs()
        G = CompleteLayerGraph(inputs[i])
        G.populate_graph()

    # FOR TREE STUFF AND TESTING
    # tree = PseudoTreeInput()
    # instance = PseudoTree(tree)
    # print(instance.G.nodes, instance.G.edges)
    # print(len(instance.G.nodes))
