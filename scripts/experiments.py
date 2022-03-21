import os
import pandas as pd

def single_run():
    sp_call = "./bin/sp set1_09-17-21 set1_09-17-21_26_0.4 26 5 1 6 500 30 80 10 1"
    os.system(sp_call)

def enum_vs_mip(setname):
    pass

def main():
    single_run()

if __name__ == "__main__":
    main()
