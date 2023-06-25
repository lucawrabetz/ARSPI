import os
import pandas as pd
import numpy as np

GOODRUNS = "~/Dropbox/aspi_good_data"
CSVNAME = "follower_testbed_final.csv"
RUNPATH = os.path.join(GOODRUNS, CSVNAME)
OUTNAME = "aspi_testbed_table-averages-enumbenders_latex.csv"
OUTPATH = os.path.join(GOODRUNS, OUTNAME)

ALLCOLUMNS = ["set_name","instance_name","nodes","arcs","k_zero","density","scenarios","budget","policies",
              # "MIP_OPTIMAL","MIP_objective","MIP_gap","MIP_time","MIP_partition",
              "BENDERS_OPTIMAL","BENDERS_objective","BENDERS_gap","BENDERS_time","BENDERS_cuts_rounds","BENDERS_partition",
              "ENUMERATION_OPTIMAL","ENUMERATION_objective","ENUMERATION_time","ENUMERATION_partition",
              # "GREEDY_objective","GREEDY_time","GREEDY_partition"
              ]

# output columns
COLUMNS = ["k_zero",
           "nodes","arcs",
           "scenarios","policies",
           # "MIP_objective","MIP_time","MIP_gap",
           "BENDERS_objective","BENDERS_time","BENDERS_gap","BENDERS_cuts_rounds",
           "ENUMERATION_objective","ENUMERATION_time",
           # "GREEDY_objective","GREEDY_time","approximation_ratio"
           ]

NUM_COLUMNS = [# "MIP_objective","MIP_time","MIP_gap",
               "BENDERS_objective","BENDERS_time","BENDERS_gap","BENDERS_cuts_rounds",
               "ENUMERATION_objective","ENUMERATION_time",
               # "GREEDY_objective","GREEDY_time","approximation_ratio"
]

# output columns for averages
AVG_COLUMNS = ["k_zero",
           "nodes","arcs",
           "scenarios","policies",
           # "MIP_time","MIP_gap",
           "BENDERS_time","BENDERS_gap","BENDERS_cuts_rounds",
           "ENUMERATION_time", "ratio",
           # "GREEDY_time","approximation_ratio"
               ]

def ms_to_s(data, index):
    data[index] = data[index].astype(float)
    data[index] = data[index].div(1000).round(2)

def round_to_dp(data, index, dp):
    data[index] = data[index].round(dp)

def round_to_int(data, index):
    data[index] = data[index].astype(int)

def convert_numerical_to_float(data):
    for col in NUM_COLUMNS:
        data[col] = data[col].astype(float)

def convert_all_times(data):
    ms_to_s(data, "BENDERS_time")
    ms_to_s(data, "ENUMERATION_time")

def round_columns(data):
    round_to_dp(data, "BENDERS_time", 2)
    round_to_dp(data, "ENUMERATION_time", 2)
    round_to_dp(data, "BENDERS_gap", 2)
    round_to_dp(data, "ratio", 2)
    round_to_int(data, "BENDERS_objective")
    round_to_int(data, "ENUMERATION_objective")

def k_new_rows_for_instance(group, max_k):
    # list of series
    new_rows = []
    grouped_by_k = group.groupby(["policies"])
    for k, df in grouped_by_k:
        if k == 0: continue
        # compute all means
        mean_df = df.mean()
        new_rows.append(mean_df)
    return new_rows

def compute_ratio_column(data):
    data["ratio"] = data["BENDERS_objective"] / data["ENUMERATION_objective"]
    return data

def take_averages(data):
    grouped_by_n = data.groupby(["nodes"])
    final_rows = []
    for name, group in grouped_by_n:
        new_rows = []
        grouped_by_followers = group.groupby(["scenarios"])
        for p, df in grouped_by_followers:
            mean_df = df.mean()
            new_rows.append(mean_df)
        final_rows.extend(new_rows)
    final_df = pd.DataFrame(final_rows)
    return final_df

if __name__=="__main__":
    run_df = pd.read_csv(RUNPATH)
    run_df.replace(",", "", regex=True, inplace=True)
    convert_numerical_to_float(run_df)
    run_df = compute_ratio_column(run_df)
    avg_df = take_averages(run_df)
    convert_all_times(avg_df)
    round_columns(avg_df)
    avg_df.to_csv(path_or_buf=OUTPATH, sep='&', columns=AVG_COLUMNS )
    # print(run_df)
