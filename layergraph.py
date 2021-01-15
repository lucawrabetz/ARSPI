import random


def layer_communicates(current):
    # add communication arcs within a layer - avoid graphs with no s-t path
    comm_arcs = []
    for i in current[:-1]:
        new_arc1 = (i, i+1)
        new_arc2 = (i+1, i)
        comm_arcs.append(new_arc1)
        comm_arcs.append(new_arc2)
    return comm_arcs


def arcs_for_current_layer(current, next, arcs_per_node):
    # add arcs for this layer - communication + outgoing
    arcs = []
    # if the next layer is t
    if len(next) == 1:
        sampled_nodes = random.sample(current, arcs_per_node)
        j = next[0]
        for i in sampled_nodes:
            new_arc = (i, j)
            arcs.append(new_arc)

    # for all other layers
    else:
        for i in current:
            sampled_nodes = random.sample(next, arcs_per_node)
            for j in sampled_nodes:
                new_arc = (i, j)
                arcs.append(new_arc)

    arcs.extend(layer_communicates(current))
    return arcs


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
    arcs_per_node = 0
    arcs = []
    n = 0
    m = 0
    t = 0

    def __init__(self, num_layerss, num_per_layerr, arcs_per_nodee):
        # only generates graph topology (no cost information in this class)
        self.num_layers = num_layerss
        self.num_per_layer = num_per_layerr
        self.arcs_per_node = arcs_per_nodee
        self.n = 2 + self.num_layers*self.num_per_layer
        self.t = self.n-1

        current_layer = [self.s]
        next_layer = [i for i in range(1, self.num_per_layer+1)]

        for i in range(self.num_layers+1):
            new_arcs = arcs_for_current_layer(
                current_layer, next_layer, self.arcs_per_node)
            self.arcs.extend(new_arcs)
            current_layer = next_layer
            if i == self.num_layers-1:
                next_layer = [self.t]
            else:
                next_layer = [(i + self.num_per_layer) for i in current_layer]

        self.m = len(self.arcs)

    def printGraph(self):
        print("n: " + str(self.n) + " - " + str([i for i in range(self.n)]))
        print("m: " + str(self.m))
        print("arcs: ")
        for i in range(self.m):
            print("     " + str(self.arcs[i]))


class TestBed:
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


# num_layers = 1
# num_per_layer = 1
# arcs_per_node = 1
# ll = 2
# samples = 2
# mu = 100
# sigma = 10

# bed = TestBed(num_layers, num_per_layer, arcs_per_node, ll, samples, mu, sigma)
# bed.G.printGraph()
# print(bed.cc)
