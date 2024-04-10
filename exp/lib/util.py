import os
import warnings
from typing import Dict, Any
import pandas as pd
from datetime import date


# GLOBAL COLUMN SET UP
# BASE COLUMN SETS (SMALL BUILDING BLOCKS)
# BASE / BUILDING BLOCKS, LISTS HARD-INITIALIZED
class feature:
    def __init__(self, name: str, default: Any):
        self.name = name
        self.default = default


COLS = {
    "name_inputs_str": [
        feature("set_name", "layer"),
        feature("instance_name", "layer_123"),
    ],
    "graph_inputs_int": [
        feature("nodes", 123),
        feature("arcs", 1231),
        feature("k_zero", 5),
    ],
    "graph_inputs_rat": [
        feature("density", 1231 / (123 * 122)),
    ],
    "cost_instance_inputs_int": [
        feature("scenarios", 5),
    ],
    "run_input_params_int": [
        feature("budget", 5),
        feature("policies", 2),
    ],
    "run_input_hyperparams_int": [
        feature("m_sym", -1),
        feature("g_sym", -1),
    ],
    "run_input_hyperparams_cat": [
        feature("subsolver", "NONE"),
    ],
    "run_input_params_cat": [
        feature("solver", "MIP"),
    ],
    "solution_outputs_rat": [
        feature("objective", 0.0),
    ],
    "solution_outputs_cat": [
        feature("unbounded", "NOT_UNBOUNDED"),
        feature("optimal", "OPTIMAL"),
        feature("partition", "0-1-1-0-0"),
    ],
    "model_outputs_int": [
        feature("cuts_rounds", 0),
        feature("cuts_added", 0),
    ],
    "model_outputs_rat": [
        feature("gap", 0.0),
    ],
    "time_outputs_rat": [
        feature("avg_cbtime", 0.0),
        feature("avg_sptime", 0.0),
        feature("time", 0.0),
    ],
    "slow_constants": [
        feature("empirical_optimal_ratio", -1),
        feature("empirical_suboptimal_ratio", -1),
        feature("best_optimal", -1),
        feature("best_objective", -1),
        feature("exact_alpha", -1),
        feature("exact_alpha_time_s", -1),
        feature("alpha_hat_one", -1),
        feature("alpha_hat_one_time_s", -1),
        feature("alpha_hat_two", -1),
        feature("alpha_hat_two_time_s", -1),
    ],
}

COLS["name_inputs"] = COLS["name_inputs_str"]
COLS["graph_inputs"] = COLS["graph_inputs_int"] + COLS["graph_inputs_rat"]
COLS["cost_instance_inputs"] = COLS["cost_instance_inputs_int"]
COLS["run_input_params"] = COLS["run_input_params_int"] + COLS["run_input_params_cat"]
COLS["run_input_hyperparams"] = (
    COLS["run_input_hyperparams_int"] + COLS["run_input_hyperparams_cat"]
)
COLS["run_inputs"] = COLS["run_input_params"] + COLS["run_input_hyperparams"]

COLS["inputs"] = (
    COLS["name_inputs"]
    + COLS["graph_inputs"]
    + COLS["cost_instance_inputs"]
    + COLS["run_inputs"]
)

COLS["model_outputs"] = COLS["model_outputs_int"] + COLS["model_outputs_rat"]

COLS["outputs_int"] = COLS["model_outputs_int"]
COLS["outputs_rat"] = COLS["model_outputs_rat"] + COLS["time_outputs_rat"]
COLS["outputs_cat"] = COLS["solution_outputs_cat"]

COLS["solution_outputs"] = COLS["solution_outputs_rat"] + COLS["solution_outputs_cat"]
COLS["outputs"] = (
    COLS["solution_outputs"] + COLS["model_outputs"] + COLS["time_outputs_rat"]
)


COLS["raw"] = COLS["inputs"] + COLS["outputs"]
# "processed" by slow computations in add_ratios.py. (persistent / db)
COLS["processed"] = COLS["raw"] + COLS["slow_constants"]
# "finished" by fast clean up in common_cleanup() (local to callstack)
COLS["time_outputs_s_rat"] = [
    feature(i.name + "_s", 0.0) for i in COLS["time_outputs_rat"]
]
COLS["logging_run_header"] = COLS["name_inputs_str"] + COLS["run_inputs"]
COLS["logging_outputs"] = COLS["solution_outputs"] + COLS["time_outputs_s_rat"]
COLS["finished"] = COLS["processed"] + COLS["time_outputs_s_rat"]

# Additional categorizations (no new columns past this point).
# Grouping - column subsets useful for common groupby operations, such as finding dups.
COLS["same_instance"] = (
    COLS["name_inputs"] + COLS["graph_inputs"] + COLS["cost_instance_inputs"]
)
# TODO REORGANIZE THESE CONCEPTS
COLS["same_run"] = COLS["same_instance"] + COLS["run_input_params_int"]
COLS["same_run_params"] = COLS["same_instance"] + COLS["run_input_params"]
COLS["same_run_hyperparams"] = COLS["inputs"]
COLS["same_run_and_outputs"] = (
    COLS["inputs"] + COLS["model_outputs"] + COLS["time_outputs_rat"]
)

# Types - column subsets partitioned by type, useful for cleaning functions - including all columns in "finished" - these subsets should be used in common_cleanup() once we are at "finished" columns.
COLS["integer"] = (
    COLS["graph_inputs_int"]
    + COLS["cost_instance_inputs_int"]
    + COLS["run_input_params_int"]
    + COLS["run_input_hyperparams_int"]
)
COLS["rational"] = (
    COLS["graph_inputs_rat"]
    + COLS["outputs_rat"]
    + COLS["time_outputs_s_rat"]
    + COLS["slow_constants"]
)

DP = {key.name: 2 for key in COLS["rational"]}

SOLVER_FLAGS = {
    "MIP": set({"m", "mip", "sp"}),
    "BENDERS": set({"b", "benders"}),
    "ENUMERATION": set({"e", "enum", "enumeration"}),
    "GREEDY": set({"g", "greedy"}),
}

SUBSOLVER_FLAGS = {
    "MIP": set({"m", "mip", "sp"}),
    "BENDERS": set({"b", "benders"}),
    "NONE": set({"n", "none", "na"}),
}

all_columns = set([])
for d in COLS.values():
    for key in d:
        all_columns.add(key)

COLLOG = {
    "pretty": {
        "set_name": "Set Name",
        "instance_name": "Instance Name",
        "nodes": "Nodes",
        "arcs": "Arcs",
        "k_zero": "Groups",
        "density": "Density",
        "scenarios": "Followers",
        "budget": "Budget",
        "policies": "Policies",
        "subsolver": "Sub Solver",
        "solver": "Solver",
        "unbounded": "Unbounded",
        "optimal": "Optimal",
        "objective": "Objective",
        "gap": "Gap",
        "time": "Running Time (ms)",
        "cuts_rounds": "Callbacks",
        "cuts_added": "Cuts Added",
        "avg_cbtime": "Average Callback Time (ms)",
        "avg_sptime": "Average Cut Separation Time (ms)",
        "partition": "Partition",
        "m_sym": "Manual Symmetry",
        "g_sym": "Gurobi Symmetry",
        "time_s": "Running Time (s)",
        "avg_cbtime_s": "Average Callback Time (s)",
        "avg_sptime_s": "Average Cut Separation Time (s)",
        "best_objective": "Best Objective",
        "best_optimal": "Best Optimal Objective",
        "empirical_suboptimal_ratio": "Empirical Suboptimal Ratio",
        "empirical_optimal_ratio": "Empirical Optimal Ratio",
        "exact_alpha": "Exact Alpha",
        "exact_alpha_time_s": "Exact Alpha Time (s)",
        "alpha_hat_one": "Alpha Hat One",
        "alpha_hat_one_time_s": "Alpha Hat One Time (s)",
        "alpha_hat_two": "Alpha Hat Two",
        "alpha_hat_two_time_s": "Alpha Hat Two Time (s)",
    },
    "compressed": {
        "set_name": "Set",
        "instance_name": "Instance",
        "nodes": "N",
        "arcs": "M",
        "k_zero": "k0",
        "density": "D",
        "scenarios": "P",
        "budget": "r0",
        "policies": "K",
        "solver": "Sol",
        "subsolver": "Sub",
        "unbounded": "Unb",
        "optimal": "Opt",
        "objective": "Obj",
        "gap": "Gap (%)",
        "time": "T (ms)",
        "cuts_rounds": "Cb",
        "cuts_added": "Cuts",
        "avg_cbtime": "CbT (ms)",
        "avg_sptime": "CutT (ms)",
        "partition": "Part",
        "m_sym": "Msym",
        "g_sym": "Gsym",
        "time_s": "T (s)",
        "avg_cbtime_s": "CbT (s)",
        "avg_sptime_s": "CutT (s)",
        "best_objective": "MaxObj",
        "best_optimal": "MaxOpt",
        "empirical_suboptimal_ratio": "rObj",
        "empirical_optimal_ratio": "rOpt",
        "exact_alpha": "a_exact",
        "exact_alpha_time_s": "a_e_T (s)",
        "alpha_hat_one": "a_hat1",
        "alpha_hat_one_time_s": "a_hat1_T (s)",
        "alpha_hat_two": "a_hat2",
        "alpha_hat_two_time_s": "a_hat2_T (s)",
    },
}

for x in all_columns:
    for names, d in COLLOG.items():
        if x not in d.keys():
            d[x] = x


def get_solver_from_flag(flag: str):
    for solver, flag_set in SOLVER_FLAGS.items():
        if flag.lower() in flag_set:
            return solver
    return None


def get_subsolver_from_flag(flag: str):
    """
    If the subsolver is "NONE", it will return the string
    "NONE", not a None value, which would indicate that
    an invalid argument was passed.
    """
    for solver, flag_set in SUBSOLVER_FLAGS.items():
        if flag.lower() in flag_set:
            return solver
    return None


def append_date(exp_name: str):
    """
    Append today's date to experiment name
    """
    today = date.today()
    date_str = today.strftime("%m_%d_%y")

    name = exp_name + "-" + date_str
    return name


def check_make_dir(path: str, i: int):
    """
    Recursively check if an experiment directory exists, or create one with the highest number
        - example - if "path" string is "/dat/experiments/test-01_29_22", and there already exist:
            - "/dat/experiments/test-01_29_22-0"
            - "/dat/experiments/test-01_29_22-1"
            - "/dat/experiments/test-01_29_22-2"
        we have to create the dir "/dat/experiments/test-01_29_22-3"
    """

    isdir = os.path.isdir(path + "-" + str(i))

    # if the directory exists, call on the next i
    if isdir:
        return check_make_dir(path, i + 1)

    # base case - create directory for given i (and return final path)
    else:
        os.mkdir(path + "-" + str(i))
        return path + "-" + str(i)


def ms_to_s(df: pd.DataFrame, col: str):
    """
    Convert column col in df from ms to s.
    """
    new_col = col + "_s"
    if df[col] is None:
        warnings.warn(
            "Input time ms col {} is None, conversion to seconds aborted.".format(col)
        )
        return
    df[new_col] = df[col].div(1000)


def add_seconds_columns(df: pd.DataFrame):
    """
    Convert all time columns to s.
    """
    for col in COLS["time_outputs_rat"]:
        ms_to_s(df, col.name)


def round_int_columns(df: pd.DataFrame):
    """
    Round all integer columns.
    """
    for col in COLS["integer"]:
        if col in df.columns:
            df[col] = df[col].astype(int)


def round_dp_columns(df: pd.DataFrame):
    """
    Round all decimal columns to DP decimal points.
    """
    for column in COLS["rational"]:
        col = column.name
        if col in df.columns:
            df[col] = df[col].astype(float)
            df[col] = df[col].round(DP[col])


def round(df: pd.DataFrame):
    round_int_columns(df)
    round_dp_columns(df)


def print_dict(d: Dict[Any, Any]):
    print("{")
    for k, v in d.items():
        print(k + ": " + v + ",")
    print("}")


def print_header(row):
    print(
        ", ".join(
            [
                COLLOG["pretty"][c.name] + ": " + str(row[c.name])
                for c in COLS["logging_run_header"]
            ]
        )
        + "."
    )


def print_finished_row(row):
    print(
        "    --> "
        + ", ".join(
            [
                f"{feature.name}: {row[feature.name]}"
                for feature in COLS["logging_outputs"]
            ]
        )
    )


# TODO: in common_cleanup, we should replace 0 objective (NOT_OPTIMAL) rows for the enumeration algorithm, to the uninterdicted shortest path value.
def common_cleanup(df):
    """
    1. Add any missing columns using defaults to meet COLS["processed"].
    2. Add quick processing columns to meet COLS["finished"].
    """
    # to processed.
    add_cols = {
        key.name: key.default
        for key in COLS["processed"]
        if key.name not in set(df.columns)
    }
    for col, val in add_cols.items():
        df[col] = val
    # to finished.
    add_seconds_columns(df)


def final_write(df, path):
    """
    Always starts with a df with "finished" columns, because common_clenup was called.
    """
    round(df)
    df.to_csv(
        path, columns=[c.name for c in COLS["processed"]], header=True, index=False
    )


def main():
    cols1 = COLS["same_run"]
    cols2 = COLS["inputs"]
    if set(cols1) == set(cols2):
        return
    else:
        print([col for col in cols1 if col not in cols2])
        print([col for col in cols2 if col not in cols1])


if __name__ == "__main__":
    main()
