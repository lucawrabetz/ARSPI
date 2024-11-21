import argparse
import copy
import os
import shutil

import networkx as nx
import numpy as np

from lib.util import check_make_dir
from lib.log_config import setup_logging
from lib.log_config import append_date

setup_logging()
import logging

# NEW_ARCS_PER_NODE per node in layer i
# e.g. now, we add between layers 1 and 2
# so the total added arcs would be N * size(layer1)
_DEFAULT_NEW_ARCS_PER_NODE: int = 1
_DEFAULT_NODES_PER_LAYER: int = 1
_DEFAULT_LAYERS: int = 1

# tree facts - splitting factor for layers 1->
_B_ELL: int = 2


def is_pseudo_tree(tree: nx.DiGraph) -> bool:
    t: int = tree.number_of_nodes() - 1
    tree_copy: nx.DiGraph = copy.deepcopy(tree)
    tree_copy.remove_node(t)
    return nx.is_tree(tree_copy)

def layer_nodes(tree: nx.DiGraph) -> dict[int, list[int]]:
    # return dict key = layer, val = list of nodes in layer
    # doing a bfs traversal with the edge list
    return dict(enumerate(nx.bfs_layers(tree, 0)))

def random_choice_wrapper(support, n, probs=None) -> list[int]:
    # bc np.random.choice returns an int if n==1
    if n == 1:
        return [int(np.random.choice(support, n, p=probs, replace=False))]
    return list(np.random.choice(support, n, p=probs, replace=False))


def new_edges_from_i(tree: nx.DiGraph, i_nodes: list[int], j_nodes: list[int], arcs_per_node: int, source_nodes: int) -> list[tuple[int]]:
    new_edges: list[int] = []

    sample: int = len(j_nodes) - _B_ELL
    prob: float = 1.0 / sample

    final_i_nodes = random_choice_wrapper(i_nodes, source_nodes)

    for i in final_i_nodes:
        existing = set(tree.adj[i])
        probs: list[int] = []
        for j in j_nodes:
            if j in existing:
                probs.append(0)
            else:
                probs.append(prob)
        new_js = random_choice_wrapper(j_nodes, arcs_per_node, probs)
        new_edges.extend([(i, j) for j in new_js])

    print("\nNEW EDGES:")
    print(new_edges)
    return new_edges


def degrade_tree(tree: nx.DiGraph, arcs_per_node: int, source_nodes: int, last_source_layer: int = 1) -> None:

    if not is_pseudo_tree(tree):
        raise ValueError("Tree is not a pseudo tree - are you sure you want to be degrading from this file?")

    layer_dict: dict[int, list[int]] = layer_nodes(tree)
    print(list(tree.edges))
    print(layer_dict)

    for ell, i_nodes in layer_dict.items():
        if ell == 0: continue
        new_edges: list[tuple[int]] = new_edges_from_i(tree, i_nodes, layer_dict[ell+1], arcs_per_node, source_nodes)
        tree.add_edges_from(new_edges)
        if ell == last_source_layer: break
    
    print("\nAFTER DEGRADING")
    print(list(tree.edges))

    return tree



def write_graph(graph: nx.DiGraph, dst_path: str, fn: str) -> None:
    return


def is_graph_file(f: str) -> bool:
    # f is a filename / basename
    return f.split(".")[-1] == "txt"


def main():
    # initially: number of degradations = absolute number of arcs added in layer 1-2
    # default number of degradations
    new_arcs_pnode: int = _DEFAULT_NEW_ARCS_PER_NODE

    # args: setname (assumed to only contain trees, will sanity check later)
    parser = argparse.ArgumentParser(
        prog="aspi_degrade_trees",
        description="Degrade pseudotrees layer by layer.",
        usage=f"Pass path to -s (required) for src (existing set of trees). Pass an int to -n (optional) for number new arcs per node per layer, default is {_DEFAULT_NEW_ARCS_PER_NODE}. Pass an integer to -l (optional) for number of layers to degrade. What happens if more layers are passed than exist? Default is {_DEFAULT_LAYERS}",
    )
    parser.add_argument(
        "-s", "--src_set", type=str, help="Path to the src set.", required=True
    )
    parser.add_argument(
        "-n", "--newarcs", type=int, help="New arcs per node in layer 1."
    )
    parser.add_argument(
        "-l", "--layers", type=int, help="Layers to degrade."
    )

    args = parser.parse_args()
    if args.newarcs:
        new_arcs_pnode = args.newarcs
    source_nodes = _DEFAULT_NODES_PER_LAYER
    src_path: str = os.path.abspath(args.src_set)
    src_sn: str = os.path.basename(src_path)
    print(src_path)
    print(src_sn)

    # create target directory target
    temp_path: str = os.path.join(os.path.dirname(src_path), f"{src_sn}deg{str(new_arcs_pnode)}")
    dst_path: str = check_make_dir(temp_path, 0, delim="_")
    print(dst_path)
    dst_sn: str = os.path.basename(dst_path)

    # for every .txt (graph) file in set
    for fn in os.listdir(src_path):
        src_fp = os.path.join(src_path, fn)
        if not os.path.isfile(src_fp):
            continue
        if not is_graph_file(fn):
            continue

        print(fn)
        dst_fn: str = fn.replace(src_sn, dst_sn)
        dst_fp: str = os.path.join(dst_path, dst_fn)
        shutil.copyfile(src_fp, dst_fp)

        if os.path.isfile(dst_fp):
            print(f"copied {fn} to {dst_fp}")

        # read the original tree
        tree = nx.read_edgelist(src_fp, nodetype=int, create_using=nx.DiGraph)

        # degrade the tree with n new arcs
        graph = degrade_tree(tree, new_arcs_pnode, source_nodes)
        write_graph(graph, dst_path, fn)

        # read all of the original cost files

        # add costs to all the cost copies for all the new arcs


if __name__ == "__main__":
    main()
