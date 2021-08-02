import layergraph
import os

# PRACTITIONER DEFINED EXPERIMENTAL INPUTS
# REMEMBER TO CHANGE SETNAME OR YOU WILL OVERWRITE THE RESULTS!!!
NUM_LAYERS = [2, 4]
DENSITY_P = [0.4, 0.7]
NUM_PER_LAYER = 3
DATANAME = "dat"
SETNAME = "set1_08-2-21"

# PATHS
LOGNAME = SETNAME + ".log"
DATAPATH = os.path.join(DATANAME, SETNAME)
LOGPATH = os.path.join(DATAPATH, LOGNAME)

if not os.path.isdir(DATAPATH):
    os.mkdir(DATAPATH)
import pdb; pdb.set_trace()

for layers in NUM_LAYERS:
    for p in DENSITY_P:
        G = layergraph.LayerGraph(layers, NUM_PER_LAYER, p)
        filename = SETNAME + "_" + str(layers) + "_" + str(p) + ".txt"
        filepath = os.path.join(DATAPATH, filename)
        # if not os.path.exists(filepath):
        #     with open(filepath, 'w'): pass
        G.checksNX(filepath)
        G.printGraph

# with open(FILEPATH, "w") as file:
#     for evaders in range(1, MAX_EVADERS+1):
#         current_results = []
#         line = ""
#         for sample in range(1, SAMPLES+1):
#             sample_results = model.M2Model(
#                 ExpBed.cc[evaders][sample], ExpBed.d, ExpBed.r_0, ExpBed.G.arcs, ExpBed.G.n)
#             current_results.append(sample_results)
#             line = line + str(sample_results[2][0]) + ","
#         file.write(line + "\n")
#         results[evaders] = current_results
