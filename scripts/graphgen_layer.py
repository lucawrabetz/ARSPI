import layergraph
import os

def main():
    # PRACTITIONER DEFINED EXPERIMENTAL INPUTS
    # REMEMBER TO CHANGE SETNAME OR YOU WILL OVERWRITE THE RESULTS!!!
    NUM_LAYERS = [3, 5, 7, 10, 12, 15, 17, 20, 50, 100, 200]
    DENSITY_P = [0.4, 0.6]
    NUM_PER_LAYER = 8
    DATNAME = "../dat"
    MODELSNAME = "../modelfiles"
    SETNAME = "set1_09-17-21"

    # PATHS
    LOGNAME = SETNAME + ".log"
    DATAPATH = os.path.join(DATNAME, SETNAME)
    MODELSPATH = os.path.join(MODELSNAME, SETNAME)
    LOGPATH = os.path.join(DATAPATH, LOGNAME)

    if not os.path.isdir(DATAPATH):
        os.mkdir(DATAPATH)
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

                    filepath = os.path.join(DATAPATH, filename)
                    has_path = G.checksNX(filepath)
                    filepath_cpp = os.path.join("dat", SETNAME, filename)

                n_list.append(n)
                filepath_list.append(filepath_cpp)
                G.printGraph()

        # instances_string = str(number_of_instances) + "\n"
        n_list_string = " ".join([str(i) for i in n_list]) + "\n"
        filepath_list_string = " ".join(filepath_list) + "\n"

        # logfile.write(instances_string)
        logfile.write(n_list_string)
        logfile.write(filepath_list_string)

if __name__ == "__main__":
    main()

