import time
import os
import pdb
import time
import sys
import argparse
import pandas as pd
import networkx as nx
from collections import defaultdict
from itertools import combinations

from lib.util import final_write, common_cleanup, COLS
from compute_alpha import (
    Row,
)

EXCLUDE_RATIO = {
    "solver": "MIP",
    "solver": "BENDERS",
    "solver": "ENUMERATION",
}
EXCLUDE_EXACT = {
    "set_name": "layer",
    "solver": "MIP",
    "solver": "BENDERS",
    "solver": "ENUMERATION",
    "policies": 0,
}
EXCLUDE_BOUNDS = {
    "set_name": "trees",
    "solver": "MIP",
    "solver": "BENDERS",
    "solver": "ENUMERATION",
    "policies": 0,
}


def get_max(df, col):
    if not df.empty:
        return df[col].max()
    else:
        return -1


def skip_row(row, exclude):
    for key, value in exclude.items():
        if row[key] == value:
            return True
    return False


def add_uninterdicted_shortest_path_column(df):
    shortest_path_objectives = {}
    uninterdicted_column = []
    for _, row in df.iterrows():
        this_instance_cols = [row[f.name] for f in COLS["same_instance"]]
        instance_key = tuple(this_instance_cols)
        if row["policies"] == 0:
            shortest_path_objectives[instance_key] = row["objective"]
    for _, row in df.iterrows():
        this_instance_cols = [row[f.name] for f in COLS["same_instance"]]
        instance_key = tuple(this_instance_cols)
        try:
            uninterdicted_column.append(shortest_path_objectives[instance_key])
        except KeyError:
            print("KeyError (no uninterdicted shortest path objective): ", instance_key)
            uninterdicted_column.append(-1)
    df["uninterdicted_shortest_path"] = uninterdicted_column
    
def add_adaptive_increment(df):
    """
    Add new column for the following ratio:
    for a row, the ratio = (objective) / (objective value of same instance with policies = 0)
    """
    df["adaptive_increment"] = df["objective"] / df["uninterdicted_shortest_path"]


def add_empirical_ratio(df):
    max_column = "objective"
    best_objective_values = []
    best_optimal_values = []
    same_run_cols = [c.name for c in COLS["same_run"]]
    for index, row in df.iterrows():
        if skip_row(row, EXCLUDE_RATIO):
            best_objective_values.append(-1)
            best_optimal_values.append(-1)
            continue
        objective_mask = df[same_run_cols].eq(row[same_run_cols]).all(axis=1)
        optimal_mask = df[same_run_cols].eq(row[same_run_cols]).all(axis=1) & (
            df["optimal"] == "OPTIMAL"
        )
        objective_df = df[objective_mask]
        optimal_df = df[optimal_mask]
        best_objective_values.append(get_max(objective_df, max_column))
        best_optimal_values.append(get_max(optimal_df, max_column))
    df["best_objective"] = best_objective_values
    df["best_optimal"] = best_optimal_values
    df["empirical_suboptimal_ratio"] = df["objective"] / df["best_objective"]
    df["empirical_optimal_ratio"] = df["objective"] / df["best_optimal"]
    df.loc[df["best_objective"] == -1, "empirical_suboptimal_ratio"] = -1
    df.loc[df["best_optimal"] == -1, "empirical_optimal_ratio"] = -1
    return df


def add_exact_alpha(df):
    print("adding exact alpha...")
    new_column = []
    time_column = []
    all_paths_total_costs_dict = {}
    for _, row in df.iterrows():
        if skip_row(row, EXCLUDE_EXACT):
            new_column.append(-1)
            time_column.append(-1)
            continue
        instance = Row(row)
        start = time.time()
        new_value = instance.compute_exact_alpha()
        duration = time.time() - start
        new_column.append(new_value)
        print("computed exact alpha value ", new_value, " in ", duration, " seconds.")
        time_column.append(duration)

    # Add the new column to the DataFrame
    df["exact_alpha"] = new_column
    df["exact_alpha_time_s"]
    return df


def add_alpha_one(df):
    print("adding alphahat 1...")
    # Iterate through each row and compute the value for the new column
    new_column = []
    times_column = []
    for _, row in df.iterrows():
        if skip_row(row, EXCLUDE_BOUNDS):
            new_column.append(-1)
            times_column.append(-1)
            continue
        instance = Row(row)
        begin_seconds = time.time()
        new_value = instance.compute_alphahat1()
        duration_seconds = time.time() - begin_seconds
        new_column.append(new_value)
        times_column.append(duration_seconds)
        print(
            "computed alphahat1 value ",
            new_value,
            " in ",
            duration_seconds,
            " seconds.",
        )
    # Add the new column to the DataFrame
    df["alpha_hat_one"] = new_column
    df["alpha_hat_one_time_s"] = times_column
    return df


def add_alpha_two(df):
    print("adding alphahat 2...")
    alphahatn_column = []
    times_column = []
    for _, row in df.iterrows():
        if skip_row(row, EXCLUDE_BOUNDS):
            alphahatn_column.append(-1)
            times_column.append(-1)
            continue
        instance = Row(row)
        begin_seconds = time.time()
        new_value = instance.compute_alphahat2()
        duration_seconds = time.time() - begin_seconds
        alphahatn_column.append(new_value)
        times_column.append(duration_seconds)
        print(
            "computed alphahat2 value ",
            new_value,
            " in ",
            duration_seconds,
            " seconds.",
        )
    df["alpha_hat_two"] = alphahatn_column
    df["alpha_hat_two_time_s"] = times_column
    return df


def add_alphas(df):
    add_exact_alpha(df)
    add_alpha_one(df)
    add_alpha_two(df)


def main():
    parser = argparse.ArgumentParser(
        description="Compute empirical ratio and alpha for aspi results csv file - original file location will be overwritten with new file."
    )
    parser.add_argument("file_path", help="Path to the CSV file")
    args = parser.parse_args()
    df = pd.read_csv(args.file_path)
    cleanup_to_processed(df)
    data_df = cleanup_to_finished(df)
    del df
    # add_empirical_ratio(data_df)
    add_uninterdicted_shortest_path_column(data_df)
    add_adaptive_increment(data_df)
    # add_alphas(data_df)
    final_write(data_df, args.file_path)


if __name__ == "__main__":
    main()
