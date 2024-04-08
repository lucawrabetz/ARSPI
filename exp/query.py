import pandas as pd
import argparse
from lib.util import *

OUTPUT_COLUMNS = [
    "instance_name",
    "policies",
    "budget",
    "optimal",
    "objective",
    "time_s",
    "solver",
    "m_sym",
    "g_sym",
    "gap",
    "exact_alpha",
]
COMPRESSED_OUTPUT_COLUMNS = [
    "set_name",
    "nodes",
    "k_zero",
    "scenarios",
    "policies",
    "budget",
    "objective",
    "time_s",
    "solver",
    "m_sym",
    "g_sym",
    "gap",
    "exact_alpha",
]
AVG_NOT_MATCH = [
    "time",
    "time_s",
    "instance_name",
    "optimal",
    "unbounded",
    "cuts_rounds",
    "cuts_added",
    "avg_cbtime",
    "avg_sptime",
    "avg_cbtime_s",
    "avg_sptime_s",
    "partition",
    "gap",
    "objective",
    "best_optimal",
    "best_objective",
    "empirical_optimal_ratio",
    "empirical_suboptimal_ratio",
]
X = 20


def filter_dataframe(df, args, solver):
    """
    Filter dataframe based on provided kwargs.
    """
    mask = pd.Series(True, index=df.index)
    for key, value in args.items():
        if (
            key != "file_path"
            and key != "average"
            and key != "solver"
            and key != "verbose"
            and value is not None
        ):
            mask = mask & (df[key] == value)
    if solver:
        mask = mask & (df["solver"] == solver)
    return df[mask]


def main():
    # TODO: (LW) this function should be refactored out, should always accept all columns from the data model.
    # If a column is added to the data model, the parser for it should automatically accept it.
    parser = argparse.ArgumentParser(description="Filter CSV file based on criteria.")
    parser.add_argument("file_path", help="Path to the CSV file")
    parser.add_argument("--set_name", type=str, help="Filter by set name")
    parser.add_argument("--nodes", type=int, help="Filter by number of nodes")
    parser.add_argument("--arcs", type=int, help="Filter by number of arcs")
    parser.add_argument("--k_zero", type=int, help="Filter by k zero")
    parser.add_argument("--scenarios", type=int, help="Filter by number of scenarios")
    parser.add_argument("--budget", type=int, help="Filter by budget")
    parser.add_argument("--policies", type=int, help="Filter by number of policies")
    parser.add_argument("--objective", type=float, help="Filter by number of policies")
    parser.add_argument("--solver", type=str, help="Filter by solver")
    parser.add_argument("--optimal", type=str, help="Filter by optimality")
    parser.add_argument("--subsolver", type=str, help="Filter by solver")
    parser.add_argument("--m_sym", type=int, help="Filter by manual symmetry")
    parser.add_argument("--g_sym", type=int, help="Filter by gurobi symmetry")
    parser.add_argument(
        "--average", action="store_true", help="Flag to average the matching rows"
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Verbose column header output"
    )

    args = parser.parse_args()
    if args.solver:
        solver = get_solver_from_flag(args.solver)
    else:
        solver = None
    if args.subsolver:
        subsolver = get_subsolver_from_flag(args.subsolver)
    else:
        subsolver = None

    if solver in ["MIP", "BENDERS"]:
        OUTPUT_COLUMNS.append("m_sym")
        OUTPUT_COLUMNS.append("g_sym")
        COMPRESSED_OUTPUT_COLUMNS.append("m_sym")
        COMPRESSED_OUTPUT_COLUMNS.append("g_sym")
    if solver == "GREEDY":
        OUTPUT_COLUMNS.append("empirical_suboptimal_ratio")
        OUTPUT_COLUMNS.append("empirical_optimal_ratio")
        COMPRESSED_OUTPUT_COLUMNS.append("empirical_suboptimal_ratio")
        COMPRESSED_OUTPUT_COLUMNS.append("empirical_optimal_ratio")
        if subsolver in ["MIP", "BENDERS"]:
            OUTPUT_COLUMNS.append("subsolver")
            COMPRESSED_OUTPUT_COLUMNS.append("subsolver")

    if solver == "BENDERS":
        OUTPUT_COLUMNS.append("cuts_rounds")
        OUTPUT_COLUMNS.append("cuts_added")
        OUTPUT_COLUMNS.append("avg_cbtime")
        OUTPUT_COLUMNS.append("avg_sptime")
        COMPRESSED_OUTPUT_COLUMNS.append("cuts_rounds")
        COMPRESSED_OUTPUT_COLUMNS.append("cuts_added")
        COMPRESSED_OUTPUT_COLUMNS.append("avg_cbtime")
        COMPRESSED_OUTPUT_COLUMNS.append("avg_sptime")

    try:
        df = pd.read_csv(args.file_path)
    except FileNotFoundError:
        print("Error: File not found.")
        return

    # COMMON CLEANUP
    common_cleanup(df)

    print("\n")
    if solver:
        print("Filtering on solver: " + solver)
    print("Arguments:\n")
    print(
        "".join(
            [
                f"    {arg}: {value}\n"
                for arg, value in vars(args).items()
                if value and arg != "average" and arg != "solver"
            ]
        )
    )

    if args.verbose:
        colname_map = COLLOG["pretty"]
    else:
        colname_map = COLLOG["compressed"]

    filtered_df = filter_dataframe(df, vars(args), solver)

    if filtered_df.empty:
        print("No matching rows found.")
    else:
        if len(filtered_df) > X:
            print(f"The filtered DataFrame has more than {X} rows.")
            choice = (
                input(
                    "Do you want to see all rows (type 'all'), just the first "
                    + str(X)
                    + " rows (type 'first'), or skip printing (type 'skip')? "
                )
                .strip()
                .lower()
            )
            print("\n")
            if choice == "all":
                print(filtered_df[OUTPUT_COLUMNS].rename(columns=colname_map))
            elif choice == "first":
                print(filtered_df.head(X)[OUTPUT_COLUMNS].rename(columns=colname_map))
            else:
                print("Skipping printing.")
        else:
            print("\n")
            print(filtered_df[OUTPUT_COLUMNS].rename(columns=colname_map))

        if args.average:
            group_cols = [
                col for col in filtered_df.columns if col not in AVG_NOT_MATCH
            ]
            compressed_df = filtered_df.groupby(group_cols).mean().reset_index()
            print("\n")
            print("Averaged: ")
            print(compressed_df[COMPRESSED_OUTPUT_COLUMNS].rename(columns=colname_map))


if __name__ == "__main__":
    main()
