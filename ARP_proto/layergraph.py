import random
import math
import networkx as nx

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


def cost_generation(num_evaders, m, mu, sigma):
    # defines a [discrete] uncertainty set of costs for some number of policies/evaders
    # samples from a normal dist.
    costs = []
    for i in range(num_evaders):
        this_evader = []
        for j in range(m):
            this_evader.append(int(abs(random.gauss(mu, sigma))))
        costs.append(this_evader)
    return costs

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


class TestBed:
    '''
    - adding costs to a layerGraph to create a full testbed
    - as of right now, not being used as I do this directly in the cpp code
    '''
    G = None
    l = 0
    samples = 0
    mu = 0
    sigma = 0
    cc = {}
    d = 0
    r_0 = 0

    def __init__(self, num_layerss, num_per_layerr, arcs_per_nodee, ll, sampless, muu, sigmaa, r_00):
        self.G = LayerGraph(num_layerss, num_per_layerr, arcs_per_nodee)
        self.l = ll
        self.samples = sampless
        self.mu = muu
        self.sigma = sigmaa
        self.d = self.sigma*2
        self.r_0 = r_00

        for evaders in range(1, self.l+1):
            current_evader_num = {}
            for i in range(1, self.samples+1):
                current_evader_num[i] = cost_generation(
                    evaders, self.G.m, self.mu, self.sigma)
            self.cc[evaders] = current_evader_num

    def writeBed(self, filename):
        with open(filename, "w") as file:
            file.write("n: " + str(self.G.n) + "\n")
            file.write("m: " + str(self.G.m) + "\n")
            file.write("R0: " + str(self.r_0) + "\n")
            file.write("max evaders: " + str(self.l) + "\n")
            file.write("samples: " + str(self.samples) + "\n")
            file.write("mu: " + str(self.mu) + "\n")
            file.write("sigma: " + str(self.sigma) + "\n")


if __name__ == "__main__":

    num_layers = 4
    num_per_layer = 8
    p = 0.7
    ll = 2
    samples = 2
    mu = 100
    sigma = 10
    r_0 = 1

    lG = LayerGraph(num_layers, num_per_layer, p)
    lG.printGraph()
    lG.checksNX('graph1.graph')
