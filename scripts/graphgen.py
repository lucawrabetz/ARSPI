import graphclasses
import os

def main():
    # PRACTITIONER DEFINED EXPERIMENTAL INPUTS
    # REMEMBER TO CHANGE SETNAME OR YOU WILL OVERWRITE THE RESULTS!!!
    N_VALUES = [10, 15, 20, 50, 100, 200]
    PR_VALUES = [0.4, 0.6]
    DATNAME = "dat"
    MODELSNAME = "modelfiles"
    SETNAME = "set3_03-24-22"

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
    number_of_instances = len(N_VALUES) * len(PR_VALUES)
    filepath_list = [] # append through the loop to ensure it matches the order
    n_list = [] # append through the loop to ensure it matches the order

    with open (LOGPATH, "w") as logfile:
        for n in N_VALUES:
            for pr in PR_VALUES:
                G = graphclasses.ErdosRenyi(n, pr)
                filename = SETNAME + "_" + str(n) + "_" + str(pr) + ".txt"
                filepath = os.path.join(DATAPATH, filename)
                G.checksNX(filepath)
                n_list.append(n)
                filepath_list.append(filepath)

        # instances_string = str(number_of_instances) + "\n"
        n_list_string = " ".join([str(i) for i in n_list]) + "\n"
        filepath_list_string = " ".join(filepath_list) + "\n"

        # logfile.write(instances_string)
        logfile.write(n_list_string)
        logfile.write(filepath_list_string)

if __name__ == "__main__":
    main()

