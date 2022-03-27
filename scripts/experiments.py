import os
import sys
import pdb
import subprocess
import math
import pandas as pd

from datetime import date

# example call:
# ./bin/sp set1_09-17-21 set1_09-17-21_26_0.4 26 5 1 6 500 30 80 10 1
# ./bin/<algorithm executable> <instance setname> <instance fullname> <n> <p> <k> <r_0> <M> <lb> <ub> <d> <costs>
# output: <objective> <runtime(ms)>

# Global namespace 
DAT = "dat"
MODELS = "modelfiles"
M = 500
LB_COST = 30 # cost distribution lower bound
UB_COST = 80 # cost distribution upper bound
INT_COST = 10 # interdiction cost increase
COLUMNS = ["setname", "name", "n", "m", "pr", "density", "k", "p", "r_0", "mip_time", "mip_obj", "mip_gap", "enum_time", "enum_obj"]
# for any given p (number of scenarios/followers) that you'd like to experiment with
# we have a set list of ks to run, usually 5 values of k (unless p is very small)
K_LISTS = {
    1: [1],
    2: [1, 2],
    3: [1, 2, 3],
    4: [1, 2, 3, 4],
    5: [1, 2, 3, 4, 5],
    6: [1, 2, 3, 4, 6],
    7: [1, 2, 4, 6, 7],
    8: [1, 2, 4, 6, 8],
    9: [1, 2, 4, 6, 9],
    10: [1, 3, 5, 8, 10],
    11: [1, 3, 5, 8, 11],
    12: [1, 3, 6, 9, 12],
    13: [1, 3, 6, 9, 13],
    14: [1, 3, 6, 9, 14],
    15: [1, 5, 8, 11, 15],
    16: [1, 5, 8, 11, 16],
    17: [1, 5, 9, 13, 17],
    18: [1, 5, 9, 13, 18],
    19: [1, 5, 9, 13, 19],
    20: [1, 5, 10, 15, 20]
}

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


def single_run(setname, name, n, p, k, r_0):
    """
    Perform a run of MIP set partitioning formulation, and enumeration algorithm on instance
    Return objective and time for both
    """
    # check for cost file to see if it needs to be generated
    # if it does, it will be generated in sp_call
    costs_filename = name + "-costs_" + str(p) + ".csv"
    costs_filepath = os.path.join(DAT, setname, costs_filename)
    costs = int(not os.path.exists(costs_filepath))

    sp_call = ["./bin/sp", setname, name, str(n), str(p), str(k), str(r_0), str(M), str(LB_COST), str(UB_COST), str(INT_COST), str(costs)]
    # no costs in enum_call because definitely generated in sp_call
    enum_call = ["./bin/enum", setname, name, str(n), str(p), str(k), str(r_0), str(M), str(LB_COST), str(UB_COST), str(INT_COST), '0']

    # .run(<executable>, <saves stdout>).<more options to encode stdout>.<throw out gurobi prints>
    sp_result = subprocess.run(sp_call, stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')[-2]
    mip_obj = float(sp_result.split(' ')[0]) * (-1) # fix negative
    mip_time = float(sp_result.split(' ')[1]) * (0.001) # convert to seconds
    mip_gap = float(sp_result.split(' ')[2])

    enum_result = subprocess.run(enum_call, stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')[-2]
    enum_obj = float(enum_result.split(' ')[0])
    enum_time = float(enum_result.split(' ')[1]) * (0.001) # convert to seconds

    full_result = [mip_time, mip_obj, mip_gap, enum_time, enum_obj]
    print(full_result)

    return full_result


def enum_vs_mip(setname, p_list):
    logfile = setname + ".log"
    logpath = os.path.join(DAT, setname, logfile)

    results = []

    with open (logpath, "r") as file:
        lines = file.readlines()
        n_list = lines[0].split(" ")
        m_list = lines[1].split(" ")
        pr_list = lines[2].split(" ")
        density_list = lines[3].split(" ")
        name_list = lines[4].split(" ")

    for i in range(len(n_list)):
        for p in p_list:
            for k in K_LISTS[p]:
                n = int(n_list[i].strip())
                m = int(m_list[i].strip())
                pr = float(pr_list[i].strip())
                density = float(density_list[i].strip())
                name = name_list[i].strip()
                r_0 = math.floor(n/4)

                single_result = [setname, name, n, m, pr, density, k, p, r_0]
                run_result = single_run(setname, name, n, p, k, r_0)
                single_result.extend(run_result)
                results.append(single_result)

    results_df = pd.DataFrame(results, columns=COLUMNS)

    temp_path = os.path.join(DAT, setname, append_date("results"))
    results_path = check_make_dir(temp_path, 0)

    csv_filename = "results" + ".csv"
    csv_path = os.path.join(results_path, csv_filename)
    results_df.to_csv(csv_path)
    print(results_df)


def main():
    """
    Inputs:
        - cpp calls require set name, full name (no paths and no file extension)
            - set name passed as arg to this script
            - full names are in the log file
    """
    # SET = sys.argv[1]
    setlist = ["first-03_27_22-0", "first-03_27_22-1", "first-03_27_22-2"]
    p_list1 = [2, 3]
    p_list2 = [4, 5]

    for SET in setlist:
        enum_vs_mip(SET, p_list1)
    for SET in setlist:
        enum_vs_mip(SET, p_list2)

if __name__ == "__main__":
    main()
