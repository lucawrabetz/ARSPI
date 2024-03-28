import os
import pdb
import time
import sys
import pandas as pd
import networkx as nx
from collections import defaultdict
from itertools import combinations
import argparse

from lib.util import *


def get_max(df, col):
    if not df.empty:
        return df[col].max()
    else:
        return -1


def add_empirical_ratio(df):
    max_column = "objective"
    best_objective_values = []
    best_optimal_values = []
    same_run_cols = [c.name for c in COLS["same_run"]]
    for index, row in df.iterrows():
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


def main():
    parser = argparse.ArgumentParser(
        description="Compute empirical ratio and alpha for aspi results csv file - original file location will be overwritten with new file."
    )
    parser.add_argument("file_path", help="Path to the CSV file")
    args = parser.parse_args()
    df = pd.read_csv(args.file_path)
    common_cleanup(df)
    add_empirical_ratio(df)
    final_write(df, args.file_path)


if __name__ == "__main__":
    main()
