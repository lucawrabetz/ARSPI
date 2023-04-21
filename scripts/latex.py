import os
import pandas as pd

GOODRUNS = "/home/luw28/good_runs"
CSVNAME = "aspi_testbed_final_run.csv"
RUNPATH = os.path.join(GOODRUNS, CSVNAME)
OUTNAME = "aspi_testbed_final_run_latex.csv"
OUTPATH = os.path.join(GOODRUNS, OUTNAME)

ALLCOLUMNS = ["set_name","instance_name","nodes","arcs","k_zero","density","scenarios","budget","policies",
              "MIP_OPTIMAL","MIP_objective","MIP_gap","MIP_time","MIP_partition",
              "BENDERS_OPTIMAL","BENDERS_objective","BENDERS_gap","BENDERS_time","BENDERS_cuts_rounds","BENDERS_partition",
              "ENUMERATION_OPTIMAL","ENUMERATION_objective","ENUMERATION_time","ENUMERATION_partition",
              "GREEDY_objective","GREEDY_time","GREEDY_partition"
              ]

COLUMNS = ["k_zero",
           "nodes","arcs",
           "scenarios","policies",
           "MIP_objective","MIP_time","MIP_gap",
           "BENDERS_objective","BENDERS_time","BENDERS_gap","BENDERS_cuts_rounds",
           "ENUMERATION_objective","ENUMERATION_time",
           "GREEDY_objective","GREEDY_time"
              ]

def ms_to_s(data, index):
    data[index] = data[index].astype(float)
    data[index] = data[index].div(1000).round(2)


def round_to_dp(data, index, dp):
    data[index] = data[index].astype(float)
    data[index] = data[index].round(dp)

if __name__=="__main__":
    run_df = pd.read_csv(RUNPATH)
    run_df.replace(",", "", regex=True, inplace=True)
    ms_to_s(run_df, "MIP_time")
    ms_to_s(run_df, "BENDERS_time")
    ms_to_s(run_df, "ENUMERATION_time")
    ms_to_s(run_df, "GREEDY_time")
    round_to_dp(run_df, "MIP_objective", 0)
    round_to_dp(run_df, "BENDERS_objective", 0)
    round_to_dp(run_df, "ENUMERATION_objective", 0)
    round_to_dp(run_df, "GREEDY_objective", 0)
    round_to_dp(run_df, "MIP_gap", 2)
    round_to_dp(run_df, "BENDERS_gap", 2)
    print(run_df)
    run_df.to_csv(path_or_buf=OUTPATH, sep="&", columns=COLUMNS)
