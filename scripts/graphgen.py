import os
import sys
import graphclasses
from datetime import date

# Globale namespace
DAT = "dat"
MODELS = "modelfiles"

def append_date(basename):
    """
    Append today's date to experiment name
    """
    today = date.today()
    date_str = today.strftime("%m_%d_%y")

    name = basename + "-" + date_str
    return name


def check_make_dir(path, i):
    """
    Recursively check if an experiment directory exists, or create one with the highest number
        - example - if "path" string is "/dat/experiments/test-01_29_22", and there already exist:
            - "/dat/experiments/test-01_29_22-0"
            - "/dat/experiments/test-01_29_22-1"
            - "/dat/experiments/test-01_29_22-2"
        we have to create the dir "/dat/experiments/test-01_29_22-3"
    """

    isdir = os.path.isdir(path + "-" + str(i))

    # if the directory exists, call on the next i
    if isdir:
        return check_make_dir(path, i + 1)

    # base case - create directory for given i (and return final path)
    else:
        os.mkdir(path + "-" + str(i))
        return path + "-" + str(i)


def main():
    """
        - Pass one arg - the desired "basename" (e.g. "luca_graph")
            - do not use the char '-' use different delimeter such as '_'
    """
    # PRACTITIONER DEFINED EXPERIMENTAL INPUTS
    N_VALUES = [10, 15, 20, 50, 100, 200]
    PR_VALUES = [0.4, 0.6]

    BASENAME = sys.argv[1]
    setname = append_date(BASENAME)
    dat_temppath = os.path.join(DAT, setname)
    models_temppath = os.path.join(MODELS, setname)

    DATPATH = check_make_dir(dat_temppath, 0)
    MODELSPATH = check_make_dir(models_temppath, 0)
    SETNAME = DATPATH.split("/")[-1]

    # PATHS
    LOGNAME = SETNAME + ".log"
    LOGPATH = os.path.join(DATPATH, LOGNAME)

    fullname_list = [] # append through the loop to ensure it matches the order
    n_list = [] # append through the loop to ensure it matches the order

    with open (LOGPATH, "w") as logfile:
        for n in N_VALUES:
            for pr in PR_VALUES:
                G = graphclasses.ErdosRenyi(n, pr)
                fullname = SETNAME + "_" + str(n) + "_" + str(pr)
                filename = fullname + ".txt"
                filepath = os.path.join(DATPATH, filename)
                G.checksNX(filepath)
                n_list.append(n)
                fullname_list.append(fullname)

        # instances_string = str(number_of_instances) + "\n"
        n_list_string = " ".join([str(i) for i in n_list]) + "\n"
        fullname_list_string = " ".join(fullname_list) + "\n"

        # logfile.write(instances_string)
        logfile.write(n_list_string)
        logfile.write(fullname_list_string)

if __name__ == "__main__":
    main()

