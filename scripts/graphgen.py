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
        - Args
            - (1) the desired "basename" (e.g. "luca_graph")
            - (2) max n
            - (3) max pr
            - do not use the char '-' use different delimeter such as '_'
    """
    N_VALUES = [10, 15, 20, 50, 100, 200, 400, 600, 800, 1000, 2000, 5000, 10000, 100000, 500000, 1000000]
    PR_VALUES = [0.3, 0.4, 0.5, 0.6, 0.7]

    # will set some default args if you don't pass them right
    if len(sys.argv) == 4:
        MAX_PR = float(sys.argv[3])
        MAX_N = int(sys.argv[2])
        BASENAME = sys.argv[1]
    elif len(sys.argv) == 3:
        MAX_PR = 0.5
        MAX_N = int(sys.argv[2])
        BASENAME = sys.argv[1]
    elif len(sys.argv) == 2:
        MAX_PR = 0.5
        MAX_N = 1000
        BASENAME = sys.argv[1]
    else:
        MAX_PR = 0.5
        MAX_N = 1000
        BASENAME = "graphs"


    setname = append_date(BASENAME)
    dat_temppath = os.path.join(DAT, setname)
    models_temppath = os.path.join(MODELS, setname)

    DATPATH = check_make_dir(dat_temppath, 0)
    MODELSPATH = check_make_dir(models_temppath, 0)
    SETNAME = DATPATH.split("/")[-1]

    # PATHS
    LOGNAME = SETNAME + ".log"
    LOGPATH = os.path.join(DATPATH, LOGNAME)

    # append through the loop to ensure it matches the order
    fullname_list = []
    n_list = []
    m_list = []
    pr_list = []
    density_list = []

    with open (LOGPATH, "w") as logfile:
        for n in N_VALUES:
            if n > MAX_N: break
            for pr in PR_VALUES:
                if pr > MAX_PR: break
                G = graphclasses.ErdosRenyi(n, pr)
                fullname = SETNAME + "-" + str(n) + "_" + str(int(10*pr))
                filename = fullname + ".txt"
                filepath = os.path.join(DATPATH, filename)
                G.checksNX(filepath)

                n_list.append(n)
                m_list.append(G.m)
                pr_list.append(pr)
                density_list.append(G.density)
                fullname_list.append(fullname)

        # instances_string = str(number_of_instances) + "\n"
        n_list_string = " ".join([str(i) for i in n_list]) + "\n"
        m_list_string = " ".join([str(i) for i in m_list]) + "\n"
        pr_list_string = " ".join([str(i) for i in pr_list]) + "\n"
        density_list_string = " ".join([str(i) for i in density_list]) + "\n"
        fullname_list_string = " ".join(fullname_list) + "\n"

        # logfile.write(instances_string)
        logfile.write(n_list_string)
        logfile.write(m_list_string)
        logfile.write(pr_list_string)
        logfile.write(density_list_string)
        logfile.write(fullname_list_string)

if __name__ == "__main__":
    main()

