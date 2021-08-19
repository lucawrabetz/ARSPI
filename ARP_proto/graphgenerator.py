import layergraph
import os

# PRACTITIONER DEFINED EXPERIMENTAL INPUTS
# REMEMBER TO CHANGE SETNAME OR YOU WILL OVERWRITE THE RESULTS!!!
NUM_LAYERS = [3, 5]
DENSITY_P = [0.5, 0.8]
NUM_PER_LAYER = 3
DATANAME = "dat"
SETNAME = "set1_08-2-21"

# PATHS
LOGNAME = SETNAME + ".log"
DATAPATH = os.path.join(DATANAME, SETNAME)
LOGPATH = os.path.join(DATAPATH, LOGNAME)

if not os.path.isdir(DATAPATH):
    os.mkdir(DATAPATH)
# import pdb; pdb.set_trace()

for layers in NUM_LAYERS:
    for p in DENSITY_P:
        has_path = False

        while not has_path:
            G = layergraph.LayerGraph(layers, NUM_PER_LAYER, p)
            filename = SETNAME + "_" + str(layers) + "_" + str(p) + ".txt"
            filepath = os.path.join(DATAPATH, filename)
            # if not os.path.exists(filepath):
            #     with open(filepath, 'w'): pass
            has_path = G.checksNX(filepath)

        G.printGraph()

# SOME OLD CODE TO GENERATE A WHOLE TESTBED (WITH COSTS, NOW IN CPP)
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
