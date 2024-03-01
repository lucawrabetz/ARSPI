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


def add_empirical_ratio(df):
    condition_columns = [col for col in SAME_RUN_COLUMNS if col != "solver"]
    max_column = "objective"

    best_objective_values = []
    best_optimal_values = []
    for index, row in df.iterrows():
        condition = df[condition_columns].eq(row[condition_columns]).all(axis=1)
        optimal_condition = df[condition_columns].eq(row[condition_columns]).all(axis=1) & (df["optimal"] == "OPTIMAL")
        filtered_df = df[condition]
        optimal_df = df[optimal_condition]
        
        if not filtered_df.empty:
            max_value = filtered_df[max_column].max().max()  # Take max value over all cells
            best_objective_values.append(max_value)
        else:
            best_objective_values.append(-1)
        if not optimal_df.empty:
            optimal_max_value = optimal_df[max_column].max().max()  # Take max value over all cells
            best_optimal_values.append(optimal_max_value)
        else:
            best_optimal_values.append(-1)


    df['best_objective'] = best_objective_values
    df['best_optimal'] = best_optimal_values
    df['empirical_suboptimal_ratio'] = df["objective"] / df["best_objective"]
    df['empirical_optimal_ratio'] = df["objective"] / df["best_optimal"]
    df.loc[df["best_objective"] == -1, 'empirical_suboptimal_ratio'] = -1
    df.loc[df["best_optimal"] == -1, 'empirical_optimal_ratio'] = -1
    return df


def main():
    parser = argparse.ArgumentParser(description='Compute empirical ratio and alpha for aspi results csv file - original file location will be overwritten with new file.')
    parser.add_argument('file_path', help='Path to the CSV file')
    args = parser.parse_args()
    df = pd.read_csv(args.file_path)
    add_empirical_ratio(df)
    import pdb; pdb.set_trace()
    df.to_csv(args.file_path, columns=FINAL_COLUMNS, index=False)


if __name__ == '__main__':
    main()
