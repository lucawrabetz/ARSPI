import pandas as pd
from lib.feature import *
from lib.util import *
from lib.paths import FINALCSVPATH

# TODO: add the capability to save the query and reload it by name or id or something.

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

# TODO: Define a class who's single job it is to receive filtered_df, the args solver and verbose, construct the temporary information about how to format printing (such as what columns to print based on what solvers etc), and then have a driver function to log the query.


def main():
    parser = FeatureArgParser("Query CSV file while filtering on features.")
    parser.add_feature_args(COLS["processed"])
    parser.add_storetruearg("--verbose")
    parser.add_storetruearg("--average")
    args = parser.parse_args()

    # print a summary of the arguments
    print("\n")
    if args.solver:
        print("Filtering on solver: " + args.solver.name)
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


    if args.solver in ["MIP", "BENDERS"]:
        OUTPUT_COLUMNS.append("m_sym")
        OUTPUT_COLUMNS.append("g_sym")
        COMPRESSED_OUTPUT_COLUMNS.append("m_sym")
        COMPRESSED_OUTPUT_COLUMNS.append("g_sym")
    if args.solver == "GREEDY":
        OUTPUT_COLUMNS.append("empirical_suboptimal_ratio")
        OUTPUT_COLUMNS.append("empirical_optimal_ratio")
        COMPRESSED_OUTPUT_COLUMNS.append("empirical_suboptimal_ratio")
        COMPRESSED_OUTPUT_COLUMNS.append("empirical_optimal_ratio")
        if subsolver in ["MIP", "BENDERS"]:
            OUTPUT_COLUMNS.append("subsolver")
            COMPRESSED_OUTPUT_COLUMNS.append("subsolver")

    if args.solver == "BENDERS":
        OUTPUT_COLUMNS.append("cuts_rounds")
        OUTPUT_COLUMNS.append("cuts_added")
        OUTPUT_COLUMNS.append("avg_cbtime")
        OUTPUT_COLUMNS.append("avg_sptime")
        COMPRESSED_OUTPUT_COLUMNS.append("cuts_rounds")
        COMPRESSED_OUTPUT_COLUMNS.append("cuts_added")
        COMPRESSED_OUTPUT_COLUMNS.append("avg_cbtime")
        COMPRESSED_OUTPUT_COLUMNS.append("avg_sptime")

    try:
        # TODO: there should be a class that reads data
        df = pd.read_csv(FINALCSVPATH)
    except FileNotFoundError:
        print("Error: File not found.")
        return

    cleanup_to_processed(df)
    data_df = cleanup_to_finished(df)
    pretty_df = cleanup_to_pretty(data_df)
    del df
    del data_df

    colname_map = dict()
    # TODO: is this dict a codesmell?
    # Construct a dictionary [Feature: string] for how they need to be printed
    for f in FEATURES:
        if args.verbose: colname_map[f.name] = f.pretty_output_name
        else: colname_map[f.name] = f.compressed_output_name

    filterer = DataFilterer(args)
    filtered_df = filterer.filter(pretty_df)

    # Printing the data:
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
