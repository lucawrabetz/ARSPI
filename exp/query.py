import pandas as pd
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
    "empirical_suboptimal_ratio",
    "exact_alpha",
    "adaptive_increment",
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
    "empirical_suboptimal_ratio",
    "exact_alpha",
    "adaptive_increment",
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
    parser = FeatureFilteringArgParser("Query CSV file while filtering on features.")
    parser.add_feature_args(COLS["processed"])
    parser.add_custom_storetruearg("--verbose")
    parser.add_custom_storetruearg("--average")

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
        df = pd.read_csv(FINALCSVPATH)
    except FileNotFoundError:
        print("Error: File not found.")
        return

    cleanup_to_processed(df)
    data_df = cleanup_to_finished(df)
    pretty_df = cleanup_to_pretty(data_df)
    del df
    del data_df

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

    filtered_df = filter_dataframe(pretty_df, vars(args), solver)

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
