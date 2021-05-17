import random

# PRACTITIONER DEFINED INPUTS
# the source node is node 0
s = 0
# define the number of layers (layer s -> layer 0; layer t -> layer num_layers + 1)
num_layers = 2
# define the number of nodes per layer
num_per_layer = 10
# define the number of outgoing arcs per node
arcs_per_node = 1
# define number of scenarios/evaders
l = 1

# COMPUTED
# compute total nodes (s + t + all layers)
n = 2 + num_layers*num_per_layer
# t is the last node
t = n-1

# helper functions


def arcs_for_current_layer(current, next, arcs_per_node=arcs_per_node):
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


# BEGIN
current_layer = [s]
next_layer = [i for i in range(1, num_per_layer+1)]
arcs = []

for i in range(num_layers+1):
    print(i, current_layer, next_layer)
    new_arcs = arcs_for_current_layer(current_layer, next_layer)
    arcs.extend(new_arcs)
    current_layer = next_layer
    if i == num_layers-1:
        next_layer = [t]
    else:
        next_layer = [(i + num_per_layer) for i in current_layer]

m = len(arcs)
print(n, m, [i for i in range(n)], arcs)
