import os
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

GOOD_DATA = "good_data_snapshots_costs"
DATA_REPO = "tests-06_06_22-0"
DAT = os.path.join(GOOD_DATA, DATA_REPO)
N = 32
P_0 = 2
K_0 = 3
P = 3
M = 75
GRAPH_NAME = DATA_REPO + "-" + str(N) + "_" + str(P_0) + "_" + str(K_0)
COST_NAME = GRAPH_NAME + "-costs_" + str(P)
GRAPH_FILE = os.path.join(DAT, GRAPH_NAME + ".txt")
COST_FILE = os.path.join(DAT, COST_NAME + ".csv")

def main():
    cost_data = pd.read_csv(COST_FILE, header=None)
    print(cost_data)
    cost_edge_labels = {}
    G = nx.read_edgelist(GRAPH_FILE)
    print(G.edges)

    with open(GRAPH_FILE) as fp:
        lines = fp.readlines()

    for a, line in enumerate(lines):
        i = str(line.split()[0])
        j = str(line.split()[1])
        edge_costs = [str(cost_data.loc[q].at[a]) for q in range(P)]
        edge_label = " ".join(edge_costs)
        edge_tuple = (i, j)
        cost_edge_labels[edge_tuple] = edge_label

    print(cost_edge_labels)
    options = {
        "pos": nx.spring_layout(G),
        "font_size": 10,
        "node_size": 250,
        "node_color": "white",
        "edgecolors": "black",
        "linewidths": 1,
        "width": 1,
        "with_labels": True,
    }


    nx.draw(G, **options)
    # nx.draw_networkx_labels(G, pos=nx.spring_layout(G))
    # nx.draw_networkx_edge_labels(G, pos=nx.spring_layout(G), edge_labels=cost_edge_labels, font_size=8)
    plt.show()


if __name__=="__main__":
    main()
