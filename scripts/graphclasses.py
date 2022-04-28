import random
import math
import networkx as nx
import pdb

def arcs_for_current_layer(current_layer, next_layer, p):
    '''
    - add arcs for this layer - outgoing
    - random arcs to nodes in all 'future layers'
    - 'next' is just all the remaining nodes in future layers
    - p is the probability that the node i is connected to any j in future layers - to simplify we use it as a proportion and convert to arcs_per_node
    '''

    new_arcs = []
    arcs_per_node = math.floor(len(next_layer) * p)

    if arcs_per_node < 1:
        arcs_per_node = 1

    for i in current_layer:
        try:
            sampled_nodes = random.sample(next_layer, arcs_per_node)
        except ValueError:
            print("SAMPLE EXCEPTION OCCURRED")
            print("len next_layer: " + str(len(next_layer)))
            print("arcs_per_node: " + str(arcs_per_node))
            print("p: " + str(p))
        for j in sampled_nodes:
            new_arc = (i, j)
            new_arcs.append(new_arc)

    return new_arcs

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
            - graphs - list of erdos-renyi graph objects
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
            if i == 0:
                temp_graphs[0] = nx.disjoint_union(sourceG, graphs[0].G)
            else:
                temp_graphs[i] = nx.disjoint_union(temp_graphs[i-1], graphs[i].G)

        pdb.set_trace()
        self.G = nx.disjoint_union(temp_graphs[-1], sinkG)

        # add the edges to connect and source and sink
        # if we want to add more edges between subgraphs we can do it here as well
        t_0 = graphs[0].t
        for i in range(len(graphs)):
            new_edges = [(self.s, n_0 * (i) + 1), (t_0 + n_0 * (i) + 1, self.t)]
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
        nx.write_edgelist(self.G, filename, data=False)


if __name__ == "__main__":
    nodes = 10
    probability = 0.1

    G1 = ErdosRenyi(nodes, probability)
    G2 = ErdosRenyi(nodes, probability)

    G3 = CompoundGraph([G1, G2])

