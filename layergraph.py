import random


def arcs_for_current_layer(current, next, arcs_per_node):
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
    return arcs


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
