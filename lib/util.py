import numbers
import pandas as pd
import numpy as np
import re
from typing import Any
from abc import ABC, abstractmethod

pandas_numeric_types = [
    pd.Int8Dtype(),
    pd.Int16Dtype(),
    pd.Int32Dtype(),
    pd.Int64Dtype(),
    pd.UInt8Dtype(),
    pd.UInt16Dtype(),
    pd.UInt32Dtype(),
    pd.UInt64Dtype(),
    pd.Float32Dtype(),
    pd.Float64Dtype(),
]
numpy_numeric_types = [
    np.int8,
    np.int16,
    np.int32,
    np.int64,
    np.uint8,
    np.uint16,
    np.uint32,
    np.uint64,
    np.float16,
    np.float32,
    np.float64,
]
ADDITIONAL_NUMERIC_TYPES = pandas_numeric_types + numpy_numeric_types


def is_numeric(dtype: type) -> bool:
    return isinstance(dtype, numbers.Number) or dtype in ADDITIONAL_NUMERIC_TYPES


def to_snakecase(s: str) -> str:
    return re.sub(r"(?<!^)(?=[A-Z])", "_", s).lower()


# GLOBAL COLUMN SET UP
# BASE FEATURE SETS (SMALL BUILDING BLOCKS)
# BASE / BUILDING BLOCKS, LISTS HARD-INITIALIZED

RAW_COLS = {
    "name_inputs_str": [
        RawFeature("set_name", str, "layer"),
        RawFeature("instance_name", str, "layer_123"),
    ],
    "graph_inputs_int": [
        RawFeature("nodes", int, 123),
        RawFeature("arcs", int, 1231),
        RawFeature("k_zero", int, 5),
    ],
    "graph_inputs_rat": [
        RawFeature("density", float, 1231 / (123 * 122)),
    ],
    "cost_instance_inputs_int": [
        RawFeature("scenarios", int, 5),
    ],
    "run_input_params_int": [
        RawFeature("budget", int, 5),
        RawFeature("policies", int, 2),
    ],
    "run_input_hyperparams_int": [
        RawFeature("m_sym", int, -1),
        RawFeature("g_sym", int, -1),
    ],
    "run_input_hyperparams_cat": [
        RawFeature("subsolver", str, "NONE"),
    ],
    "run_input_params_cat": [
        RawFeature("solver", str, "MIP"),
    ],
    "solution_outputs_rat": [
        RawFeature("objective", float, 0.0),
    ],
    "solution_outputs_cat": [
        RawFeature("unbounded", str, "NOT_UNBOUNDED"),
        RawFeature("optimal", str, "OPTIMAL"),
        RawFeature("partition", str, "0-1-1-0-0"),
    ],
    "model_outputs_int": [
        RawFeature("cuts_rounds", int, 0),
        RawFeature("cuts_added", int, 0),
    ],
    "model_outputs_rat": [
        RawFeature("gap", int, 0.0),
    ],
    "time_outputs_rat": [
        RawFeature("avg_cbtime", int, 0.0),
        RawFeature("avg_sptime", int, 0.0),
        RawFeature("time", int, 0.0),
    ],
}


class MatchConditionValue(Feature, ABC):
    MATCH: list[Feature]
    CONDITIONS: list[tuple[str, Any]]
    TIE_COL: str = "objective"
    TIE_MAX: bool = True

    @abstractmethod
    def add_column(self, df: pd.DataFrame) -> None:
        pass


class SameInstanceUninterdicted(MatchConditionValue):
    MATCH = []
    CONDITIONS = [("budget", 0), ("policies", 0), ("optimal", 0)]

    def add_column(
        self,
        df: pd.DataFrame,
    ) -> None:
        """
        Add new column df[name],
        where name = self.name,
        where df[name] = max_{r in sub} {r[TIE_COL]},
        where sub = {r in df : r[f] == df[f] for f in MATCH && r[k] == v for (k, v) in CONDITIONS}.
        If TIE_MAX = False, use min instead of max.
        """
        new_col = []
        mask_colnames = [c.name for c in self.MATCH]
        for _, row in df.iterrows():
            mask = df[mask_colnames].eq(row[mask_colnames]).all(axis=1)
            for key, val in self.CONDITIONS:
                mask = mask & (df[key] == val)
            mask_df = df[mask]
            new_col.append(get_minmax(mask_df, tie_col, tie_max))
        df[self.name] = new_col


# groups not found outside this file
COLS = RAW_COLS
COLS["slow_constants"] = list(
    [
        Feature("empirical_optimal_ratio", int, -1),
        Feature("empirical_suboptimal_ratio", int, -1),
        Feature("best_optimal", int, -1),
        Feature("best_objective", int, -1),
        Feature("same_instance_uninterdicted", int, -1),
    ]
)
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
    Feature(i.name + "_s", float, 0.0) for i in COLS["time_outputs_rat"]
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
    },
}

for x in all_columns:
    for names, d in COLLOG.items():
        if x not in d.keys():
            d[x] = x


def get_solver_from_flag(flag):
    for solver, flag_set in SOLVER_FLAGS.items():
        if flag.lower() in flag_set:
            return solver
    return None


def get_subsolver_from_flag(flag):
    """
    If the subsolver is "NONE", it will return the string
    "NONE", not a None value, which would indicate that
    an invalid argument was passed.
    """
    for solver, flag_set in SUBSOLVER_FLAGS.items():
        if flag.lower() in flag_set:
            return solver
    return None


def ms_to_s(df, col):
    """
    Convert column col in df from ms to s.
    """
    new_col = col + "_s"
    df[new_col] = df[col].div(1000)


def add_seconds_columns(df):
    """
    Convert all time columns to s.
    """
    for col in COLS["time_outputs_rat"]:
        ms_to_s(df, col.name)


def round_int_columns(df):
    """
    Round all integer columns.
    """
    for col in COLS["integer"]:
        if col in df.columns:
            df[col] = df[col].astype(int)


def round_dp_columns(df):
    """
    Round all decimal columns to DP decimal points.
    """
    for column in COLS["rational"]:
        col = column.name
        if col in df.columns:
            df[col] = df[col].astype(float)
            df[col] = df[col].round(DP[col])


def round(df):
    round_int_columns(df)
    round_dp_columns(df)


def print_dict(d):
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
                f"{Feature.name}: {row[Feature.name]}"
                for Feature in COLS["logging_outputs"]
            ]
        )
    )


def get_minmax(
    df: pd.DataFrame | pd.Series, col: str, maximum: bool = True
) -> float | Any | pd.Series:
    datatype = df[col].dtype.type
    if not df.empty and is_numeric(datatype):
        if maximum:
            return df[col].max()
        else:
            return df[col].min()
    else:
        return -1


def add_new_matchcondition_column(
    name: str,
    df: pd.DataFrame,
    match: list[Feature],
    conditions: list[tuple[str, Any]],
    tie_col: str = "objective",
    tie_max: bool = True,
) -> None:
    """
    Add new column df[name],
    where df[name] = max_{r in sub} {r[tie_col]},
    where sub = {r in df : r[f] == df[f] for f in match && r[k] == v for (k, v) in conditions}.
    If tie_max = False, use min instead of max.
    """
    new_col = []
    mask_colnames = [c.name for c in match]
    for _, row in df.iterrows():
        mask = df[mask_colnames].eq(row[mask_colnames]).all(axis=1)
        for key, val in conditions:
            mask = mask & (df[key] == val)
        mask_df = df[mask]
        new_col.append(get_minmax(mask_df, tie_col, tie_max))
    df[name] = new_col


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
