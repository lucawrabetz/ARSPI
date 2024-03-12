import os
import pandas as pd
import argparse
from lib.util import *


SET_NAME_MAPPINGS = {
    "experiment_allalgorithms": "layer",
    "experiment_bendersenum": "layer",
}

INPUT_COLUMN_MAPPING = {
    "set_name": "set_name",
    "instance_name": "instance_name",
    "nodes": "nodes",
    "arcs": "arcs",
    "k_zero": "k_zero",
    "density": "density",
    "scenarios": "scenarios",
    "budget": "budget",
    "policies": "policies",
}

MIP_DROP_COLUMNS = {
    "BENDERS_OPTIMAL": "",
    "BENDERS_objective": "",
    "BENDERS_gap": "",
    "BENDERS_time": "",
    "BENDERS_cuts_rounds": "",
    "BENDERS_partition": "",
    "ENUMERATION_OPTIMAL": "",
    "ENUMERATION_objective": "",
    "ENUMERATION_time": "",
    "ENUMERATION_partition": "",
    "GREEDY_objective": "",
    "GREEDY_time": "",
    "GREEDY_partition": "",
}

BENDERS_DROP_COLUMNS = {
    "MIP_OPTIMAL": "",
    "MIP_objective": "",
    "MIP_gap": "",
    "MIP_time": "",
    "MIP_partition": "",
    "ENUMERATION_OPTIMAL": "",
    "ENUMERATION_objective": "",
    "ENUMERATION_time": "",
    "ENUMERATION_partition": "",
    "GREEDY_objective": "",
    "GREEDY_time": "",
    "GREEDY_partition": "",
}

ENUMERATION_DROP_COLUMNS = {
    "MIP_OPTIMAL": "",
    "MIP_objective": "",
    "MIP_gap": "",
    "MIP_time": "",
    "MIP_partition": "",
    "BENDERS_OPTIMAL": "",
    "BENDERS_objective": "",
    "BENDERS_gap": "",
    "BENDERS_time": "",
    "BENDERS_cuts_rounds": "",
    "BENDERS_partition": "",
    "GREEDY_objective": "",
    "GREEDY_time": "",
    "GREEDY_partition": "",
}

GREEDY_DROP_COLUMNS = {
    "MIP_OPTIMAL": "",
    "MIP_objective": "",
    "MIP_gap": "",
    "MIP_time": "",
    "MIP_partition": "",
    "BENDERS_OPTIMAL": "",
    "BENDERS_objective": "",
    "BENDERS_gap": "",
    "BENDERS_time": "",
    "BENDERS_cuts_rounds": "",
    "BENDERS_partition": "",
    "ENUMERATION_OPTIMAL": "",
    "ENUMERATION_objective": "",
    "ENUMERATION_time": "",
    "ENUMERATION_partition": "",
}

MIP_COLUMN_MAPPING = {
    "MIP_OPTIMAL": "optimal",
    "MIP_objective": "objective",
    "MIP_gap": "gap",
    "MIP_time": "time",
    "MIP_partition": "partition",
}

MIP_COLUMNS_TO_ADD = {
    "solver": "MIP",
    "unbounded": "NOT_UNBOUNDED",
    "cuts_rounds": 0,
    "cuts_added": 0,
    "avg_cbtime": 0,
    "avg_sptime": 0,
    "m_sym": 0,
    "g_sym": -1,
}

BENDERS_COLUMN_MAPPING = {
    "BENDERS_OPTIMAL": "optimal",
    "BENDERS_objective": "objective",
    "BENDERS_gap": "gap",
    "BENDERS_time": "time",
    "BENDERS_cuts_rounds": "cuts_rounds",
    "BENDERS_partition": "partition",
}

BENDERS_COLUMNS_TO_ADD = {
    "solver": "BENDERS",
    "unbounded": "NOT_UNBOUNDED",
    "cuts_added": 0,
    "avg_cbtime": 0,
    "avg_sptime": 0,
    "m_sym": 0,
    "g_sym": -1,
}

ENUMERATION_COLUMN_MAPPING = {
    "ENUMERATION_OPTIMAL": "optimal",
    "ENUMERATION_objective": "objective",
    "ENUMERATION_time": "time",
    "ENUMERATION_partition": "partition",
}

ENUMERATION_COLUMNS_TO_ADD = {
    "solver": "ENUMERATION",
    "unbounded": "NOT_UNBOUNDED",
    "cuts_rounds": 0,
    "cuts_added": 0,
    "avg_cbtime": 0,
    "avg_sptime": 0,
    "m_sym": 0,
    "g_sym": -1,
    "gap": -1,
}

GREEDY_COLUMN_MAPPING = {
    "GREEDY_objective": "objective",
    "GREEDY_time": "time",
    "GREEDY_partition": "partition",
}

GREEDY_COLUMNS_TO_ADD = {
    "solver": "GREEDY",
    "optimal": "NOT_OPTIMAL",
    "unbounded": "NOT_UNBOUNDED",
    "cuts_rounds": 0,
    "cuts_added": 0,
    "avg_cbtime": 0,
    "avg_sptime": 0,
    "m_sym": 0,
    "g_sym": -1,
    "gap": -1,
}

SKIP_VAL = {
    "policies": 0,
}


def extract_mip(df):
    mapping = INPUT_COLUMN_MAPPING.copy()
    mapping.update(MIP_COLUMN_MAPPING)
    mapping.update(MIP_DROP_COLUMNS)
    new_columns = {
        original_name: new_name
        for original_name, new_name in mapping.items()
        if new_name != ""
    }
    new_df = df.rename(columns=new_columns)
    for key, val in MIP_COLUMNS_TO_ADD.items():
        new_df[key] = val
    return new_df


def extract_benders(df):
    mapping = INPUT_COLUMN_MAPPING.copy()
    mapping.update(BENDERS_COLUMN_MAPPING)
    mapping.update(BENDERS_DROP_COLUMNS)
    new_columns = {
        original_name: new_name
        for original_name, new_name in mapping.items()
        if new_name != ""
    }
    new_df = df.rename(columns=new_columns)
    for key, val in BENDERS_COLUMNS_TO_ADD.items():
        new_df[key] = val
    return new_df


def extract_enumeration(df):
    mapping = INPUT_COLUMN_MAPPING.copy()
    mapping.update(ENUMERATION_COLUMN_MAPPING)
    mapping.update(ENUMERATION_DROP_COLUMNS)
    new_columns = {
        original_name: new_name
        for original_name, new_name in mapping.items()
        if new_name != ""
    }
    new_df = df.rename(columns=new_columns)
    for key, val in ENUMERATION_COLUMNS_TO_ADD.items():
        new_df[key] = val
    return new_df


def extract_greedy(df):
    mapping = INPUT_COLUMN_MAPPING.copy()
    mapping.update(GREEDY_COLUMN_MAPPING)
    mapping.update(GREEDY_DROP_COLUMNS)
    new_columns = {
        original_name: new_name
        for original_name, new_name in mapping.items()
        if new_name != ""
    }
    new_df = df.rename(columns=new_columns)
    for key, val in GREEDY_COLUMNS_TO_ADD.items():
        new_df[key] = val
    return new_df


def replace_setnames(df):
    for key, val in SET_NAME_MAPPINGS.items():
        mask = df["set_name"] == key
        df.loc[mask, "set_name"] = val


def replace_instancenames(df):
    for key, val in SET_NAME_MAPPINGS.items():
        mask = df["instance_name"].str.contains(key)
        df.loc[mask, "instance_name"] = df.loc[mask, "instance_name"].str.replace(
            key, val
        )


def final_changes(df):
    df["optimal"] = df["optimal"].replace("NOT OPTIMAL", "NOT_OPTIMAL")
    for row, val in SKIP_VAL.items():
        df = df[df[row] != val]
    replace_setnames(df)
    replace_instancenames(df)
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Receive dataframe with old column format, extract single rows in new column format and write rows to a new dataframe."
    )
    parser.add_argument(
        "existing", metavar="X", type=str, nargs=1, help="existing dataframe to add to"
    )
    parser.add_argument(
        "solver", metavar="S", type=str, nargs=1, help="solver to extract"
    )
    args = parser.parse_args()
    df = pd.read_csv(args.existing[0])
    solver = args.solver[0]
    if solver == "M" or solver == "m":
        df = extract_mip(df)
    if solver == "B" or solver == "b":
        df = extract_benders(df)
    if solver == "E" or solver == "e":
        df = extract_enumeration(df)
    if solver == "G" or solver == "g":
        df = extract_greedy(df)
    df = final_changes(df)
    results_path = args.existing[0].split(".")[0] + "-split-" + solver + ".csv"
    df.to_csv(results_path, header=True, index=False, columns=COLS["processed"])


if __name__ == "__main__":
    main()
