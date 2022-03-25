import os
import sys
import math
import pandas as pd

# example call:
# ./bin/sp set1_09-17-21 set1_09-17-21_26_0.4 26 5 1 6 500 30 80 10 1
# ./bin/<algorithm executable> <instance setname> <instance fullname> <n> <p> <k> <r_0> <M> <lb> <ub> <d> <costs>
# output: <objective> <runtime(ms)>

# Global namespace 
DAT = "dat"
MODELS = "modelfiles"

def single_run(setname, name, n, p, k, costs=1):
    """
    Please pass costs as int 0 or 1, not a bool
    """
    r_0 = math.floor(n/4)
    sp_call = "./bin/sp " + setname + " " + name + " " + str(n) + " " + str(p) + " " + str(k) + " " + str(r_0) + " 500 30 80 10 " + str(costs)
    # no costs in enum_call because definitely generated in sp_call
    enum_call = "./bin/sp " + setname + " " + name + " " + str(n) + " " + str(p) + " " + str(k) + " " + str(r_0) + " 500 30 80 10 0"

    print(sp_call)
    os.system(sp_call)
    print(enum_call)
    os.system(enum_call)

def enum_vs_mip(setname):
    logfile = setname + ".log"
    logpath = os.path.join(DAT, setname, logfile)
    k = 2
    p = 5

    with open (logpath, "r") as file:
        lines = file.readlines()
        n_list = lines[0].split(" ")
        name_list = lines[1].split(" ")

    for i in range(len(n_list)):
        single_run(setname, name_list[i].strip(), int(n_list[i].strip()), p, k)


def main():
    """
    Inputs:
        - cpp calls require set name, full name (no paths and no file extension)
            - set name passed as arg to this script
            - full names are in the log file
    """
    SET = sys.argv[1]
    enum_vs_mip(SET)

if __name__ == "__main__":
    main()
