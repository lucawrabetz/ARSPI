# GLOBAL COLUMN SET UP

# BASE COLUMN SETS (SMALL BUILDING BLOCKS)
# BASE / BUILDING BLOCKS, LISTS HARD-INITIALIZED
COLS = {
    "name_inputs": [
        "set_name",
        "instance_name",
    ],
    "graph_inputs": [
        "nodes",
        "arcs",
        "k_zero",
        "density",
    ],
    "cost_instance_inputs": [
        "scenarios",
    ],
    "run_inputs": [
        "budget",
        "policies",
        "solver",
        "m_sym",
        "g_sym",
    ],
    "solution_outputs": [
        "unbounded",
        "optimal",
        "objective",
        "partition",
    ],
    "model_outputs": [
        "gap",
        "cuts_rounds",
        "cuts_added",
    ],
    "time_outputs": [
        "avg_cbtime",
        "avg_sptime",
        "time",
    ],
    "slow_constants": [
        "empirical_optimal_ratio",
        "empirical_suboptimal_ratio",
        "best_optimal",
        "best_objective",
    ],
}

COLS["inputs"] = (
    COLS["name_inputs"]
    + COLS["graph_inputs"]
    + COLS["cost_instance_inputs"]
    + COLS["run_inputs"]
)
COLS["outputs"] = (
    COLS["solution_outputs"] + COLS["model_outputs"] + COLS["time_outputs"]
)
COLS["raw"] = COLS["inputs"] + COLS["outputs"]
# "processed" by slow computations in add_ratios.py. (persistent)
COLS["processed"] = COLS["raw"] + COLS["slow_constants"]
# "finished" by fast clean up in common_cleanup() (local)
COLS["time_outputs_s"] = [i + "_s" for i in COLS["time_outputs"]]
COLS["finished"] = COLS["processed"] + COLS["time_outputs_s"]

# Additional categorizations (no new columns past this point).
additional_rational_cols = ["density"]
COLS["rational"] = (
    additional_rational_cols
    + COLS["outputs"]
    + COLS["time_outputs_s"]
    + COLS["slow_constants"]
)
COLS["same_run"] = [c for c in COLS["inputs"] if c not in ["m_sym", "g_sym"]]
COLS["same_parametrized_run"] = COLS["inputs"]
COLS["same_instance"] = [col for col in COLS["same_run"] if col != "solver"]

additional_int_columns = ["nodes", "arcs", "k_zero"]
INT_COLUMNS = additional_int_columns + COLS["cost_instance_inputs"] + COLS["run_inputs"]

DP = {key: 2 for key in COLS["rational"]}

SOLVER_FLAGS = {
    "MIP": set({"m", "mip", "sp"}),
    "BENDERS": set({"b", "benders"}),
    "ENUMERATION": set({"e", "enum", "enumeration"}),
    "GREEDY": set({"g", "greedy"}),
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
    for col in COLS["time_outputs"]:
        ms_to_s(df, col)


def round_int_columns(df):
    """
    Round all integer columns.
    """
    for col in INT_COLUMNS:
        if col in df.columns:
            df[col] = df[col].astype(int)


def round_dp_columns(df):
    """
    Round all decimal columns to DP decimal points.
    """
    for col in COLS["rational"]:
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
