import layergraph
import os

def main():
    # PRACTITIONER DEFINED EXPERIMENTAL INPUTS
    # REMEMBER TO CHANGE SETNAME OR YOU WILL OVERWRITE THE RESULTS!!!
    NUM_LAYERS = [3, 5, 10, 20, 50, 100, 200]
    DENSITY_P = [0.2]
    NUM_PER_LAYER = 4
    DATNAME = "../dat"
    MODELSNAME = "../modelfiles"
    SETNAME = "set1_09-01-21"

    # PATHS
    LOGNAME = SETNAME + ".log"
    DATAPATH = os.path.join(DATANAME, SETNAME)
    DATAPATH2 = os.path.join(DATANAME2, SETNAME)
    MODELSPATH = os.path.join(MODELSNAME, SETNAME)
    LOGPATH = os.path.join(DATAPATH, LOGNAME)

    if not os.path.isdir(DATAPATH):
        os.mkdir(DATAPATH)
    if not os.path.isdir(DATAPATH2):
        os.mkdir(DATAPATH2)
    if not os.path.isdir(MODELSPATH):
        os.mkdir(MODELSPATH)
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

                    # filepath 1 is for python, 2 is for c++
                    filepath = os.path.join(DATAPATH, filename)
                    filepath2 = os.path.join(DATAPATH2, filename)
                    has_path = G.checksNX(filepath)

                n_list.append(n)
                filepath_list.append(filepath2)
                G.printGraph()

        # instances_string = str(number_of_instances) + "\n"
        n_list_string = " ".join([str(i) for i in n_list]) + "\n"
        filepath_list_string = " ".join(filepath_list) + "\n"

        # logfile.write(instances_string)
        logfile.write(n_list_string)
        logfile.write(filepath_list_string)

if __name__ == "__main__":
    main()

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
