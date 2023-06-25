import os
import pandas as pd
import numpy as np

GOODRUNS = "~/Dropbox/aspi_good_data"
CSVNAME = "aspi_testbed_final_run.csv"
RUNPATH = os.path.join(GOODRUNS, CSVNAME)
OUTNAME = "aspi_testbed_table-averages-1_latex.csv"
OUTPATH = os.path.join(GOODRUNS, OUTNAME)

ALLCOLUMNS = ["set_name","instance_name","nodes","arcs","k_zero","density","scenarios","budget","policies",
              "MIP_OPTIMAL","MIP_objective","MIP_gap","MIP_time","MIP_partition",
              "BENDERS_OPTIMAL","BENDERS_objective","BENDERS_gap","BENDERS_time","BENDERS_cuts_rounds","BENDERS_partition",
              "ENUMERATION_OPTIMAL","ENUMERATION_objective","ENUMERATION_time","ENUMERATION_partition",
              "GREEDY_objective","GREEDY_time","GREEDY_partition"
              ]

# output columns
COLUMNS = ["k_zero",
           "nodes","arcs",
           "scenarios","policies",
           "MIP_objective","MIP_time","MIP_gap",
           "BENDERS_objective","BENDERS_time","BENDERS_gap","BENDERS_cuts_rounds",
           "ENUMERATION_objective","ENUMERATION_time",
           "GREEDY_objective","GREEDY_time","approximation_ratio"]

NUM_COLUMNS = ["MIP_objective","MIP_time","MIP_gap","MIP_adaptiveincrement",
               "BENDERS_objective","BENDERS_time","BENDERS_gap","BENDERS_cuts_rounds","BENDERS_adaptiveincrement",
               "ENUMERATION_objective","ENUMERATION_time","ENUMERATION_adaptiveincrement",
               "GREEDY_objective","GREEDY_time","approximation_ratio","GREEDY_adaptiveincrement"]

# output columns for averages
AVG_COLUMNS = ["k_zero",
           "nodes","arcs",
           "scenarios","policies",
           "MIP_time","MIP_gap","MIP_adaptiveincrement",
           "BENDERS_time","BENDERS_gap","BENDERS_cuts_rounds","BENDERS_adaptiveincrement",
           "ENUMERATION_time","ENUMERATION_adaptiveincrement",
           "GREEDY_time","approximation_ratio","GREEDY_adaptiveincrement"]

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
    ms_to_s(data, "MIP_time")
    ms_to_s(data, "BENDERS_time")
    ms_to_s(data, "ENUMERATION_time")
    ms_to_s(data, "GREEDY_time")

def round_columns(data):
    round_to_dp(data, "MIP_time", 2)
    round_to_dp(data, "BENDERS_time", 2)
    round_to_dp(data, "ENUMERATION_time", 2)
    round_to_dp(data, "GREEDY_time", 2)
    round_to_dp(data, "MIP_gap", 2)
    round_to_dp(data, "BENDERS_gap", 2)
    round_to_dp(data, "approximation_ratio", 2)
    round_to_dp(data, "MIP_adaptiveincrement", 2)
    round_to_dp(data, "BENDERS_adaptiveincrement", 2)
    round_to_dp(data, "ENUMERATION_adaptiveincrement", 2)
    round_to_dp(data, "GREEDY_adaptiveincrement", 2)
    round_to_int(data, "MIP_objective")
    round_to_int(data, "BENDERS_objective")
    round_to_int(data, "ENUMERATION_objective")
    round_to_int(data, "GREEDY_objective")

def compute_approximation_ratio(data):
    data["approximation_ratio"] = data["GREEDY_objective"] / data["ENUMERATION_objective"]
    data.loc[data["approximation_ratio"] < 0, "approximation_ratio"] = np.nan

def add_shortest_path_column(data):
    data["shortest_path"] = data[data["instance"]]

def add_shortest_path_column(df):
    # Find the MIP_objective value where policies = 0 for each instance_name
    shortest_paths = df.loc[df['policies'] == 0, ['instance_name', 'MIP_objective']]
    # Create a dictionary mapping instance_name to the corresponding shortest path
    shortest_paths_dict = shortest_paths.set_index('instance_name')['MIP_objective'].to_dict()
    # Add the 'shortest_path' column based on the matching instance_name
    df['shortest_path'] = df['instance_name'].map(shortest_paths_dict)
    return df

def compute_adaptive_increment(data):
    # Create shortest_path column
    data = add_shortest_path_column(data)
    # Add new columns
    data['MIP_adaptiveincrement'] = data['MIP_objective'] / data['shortest_path']
    data.loc[data['MIP_OPTIMAL'] == 'NOT OPTIMAL', 'MIP_adaptiveincrement'] = np.nan
    data['BENDERS_adaptiveincrement'] = data['BENDERS_objective'] / data['shortest_path']
    data.loc[data['BENDERS_OPTIMAL'] == 'NOT OPTIMAL', 'BENDERS_adaptiveincrement'] = np.nan
    data['ENUMERATION_adaptiveincrement'] = data['ENUMERATION_objective'] / data['shortest_path']
    data.loc[data['ENUMERATION_OPTIMAL'] == 'NOT OPTIMAL', 'ENUMERATION_adaptiveincrement'] = np.nan
    data['GREEDY_adaptiveincrement'] = data['GREEDY_objective'] / data['shortest_path']
    return data

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

def take_averages(data):
    # for every instance set_name make a new row with instance_name with no "id"
    grouped_by_n = data.groupby(["nodes"])
    final_rows = []
    for name, group in grouped_by_n:
        max_k = grouped_by_n["policies"].max()
        new_rows = k_new_rows_for_instance(group, max_k)
        final_rows.extend(new_rows)
    final_df = pd.DataFrame(final_rows)
    return final_df

if __name__=="__main__":
    run_df = pd.read_csv(RUNPATH)
    run_df.replace(",", "", regex=True, inplace=True)
    compute_approximation_ratio(run_df)
    run_df = compute_adaptive_increment(run_df)
    convert_numerical_to_float(run_df)
    avg_df = take_averages(run_df)
    print(avg_df.columns)
    convert_all_times(avg_df)
    round_columns(avg_df)
    avg_df.to_csv(path_or_buf=OUTPATH, sep='&', columns=AVG_COLUMNS )
    # print(run_df)
