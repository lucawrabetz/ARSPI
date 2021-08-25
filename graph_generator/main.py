import layergraph
import os

# PRACTITIONER DEFINED EXPERIMENTAL INPUTS
# REMEMBER TO CHANGE SETNAME OR YOU WILL OVERWRITE THE RESULTS!!!
NUM_LAYERS = [3, 5]
DENSITY_P = [0.5]
NUM_PER_LAYER = 3
DATANAME = "../dat"
SETNAME = "set2_08-24-21"

# PATHS
LOGNAME = SETNAME + ".log"
DATAPATH = os.path.join(DATANAME, SETNAME)
LOGPATH = os.path.join(DATAPATH, LOGNAME)

if not os.path.isdir(DATAPATH):
    os.mkdir(DATAPATH)
# import pdb; pdb.set_trace()

# OTHER DATA TO WRITE FOR CPP (WILL WRITE TO LOGFILE)
number_of_instances = len(NUM_LAYERS) * len(DENSITY_P)
filepath_list = [] # append through the loop to ensure it matches the order
n_list = [] # append through the loop to ensure it matches the order

with open (LOGPATH, "w") as logfile:
    for layers in NUM_LAYERS:
        for p in DENSITY_P:
            has_path = False
            n = (2+layers*NUM_PER_LAYER)

            while not has_path:
                G = layergraph.LayerGraph(layers, NUM_PER_LAYER, p)
                filename = SETNAME + "_" + str(n) + "_" + str(p) + ".txt"
                filepath = os.path.join(DATAPATH, filename)
                has_path = G.checksNX(filepath)

            n_list.append(n)
            filepath_list.append(filepath)
            G.printGraph()

    instances_string = str(number_of_instances) + "\n"
    n_list_string = " ".join([str(i) for i in n_list]) + "\n"
    filepath_list_string = " ".join(filepath_list) + "\n"

    logfile.write(instances_string)
    logfile.write(n_list_string)
    logfile.write(filepath_list_string)

# SOME OLD CODE TO GENERATE A WHOLE TESTBED (WITH COSTS USING THE EXPBED CLASS, NOW IN CPP)
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
