import os
import sys
import pandas as pd

# example call:
# ./bin/sp set1_09-17-21 set1_09-17-21_26_0.4 26 5 1 6 500 30 80 10 1
# ./bin/<algorithm executable> <instance setname> <instance fullname> <n> <p> <k> <> <M> <lb> <ub> <d> <costs>
# output: <objective> <runtime(ms)>

# Global namespace 
DAT = "dat"
MODELS = "modelfiles"

def single_run(setname, instancename, costs=1):
    """
    Please pass costs as int 0 or 1, not a bool
    """
    sp_call = "./bin/sp " + setname + " " + instancename + " 100 5 1 6 500 30 80 10 " + str(costs)
    print(sp_call)
    # no costs in enum_call because definitely generated in sp_call
    enum_call = "./bin/enum " + setname + " " + instancename + " 100 5 1 6 500 30 80 10 0"
    os.system(sp_call)
    os.system(enum_call)

def enum_vs_mip(setname):
    pass

def main():
    # SET = sys.argv[1]
    SET = "set1_03-24-22"
    INSTANCE = SET + "_100_0.6"
    # enum_vs_mip(SET)
    single_run(SET, INSTANCE)

if __name__ == "__main__":
    main()
