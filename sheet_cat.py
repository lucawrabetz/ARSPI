import pandas as pd
import argparse

from lib.util import *


def check_columns(A, B):
    if set(A.columns) != set(B.columns):
        print(A.columns)
        print(B.columns)
        raise ValueError(
            "Columns of DataFrame A do not match columns of DataFrame B. Cannot concatenate."
        )


def concatenate_dataframes(A, B):
    check_columns(A, B)
    combined_df = pd.concat([A, B], ignore_index=True, sort=False)
    match = [col for col in COLS["same_parametrized_run"]]
    duplicates = combined_df[combined_df.duplicated(subset=match, keep=False)]
    if not duplicates.empty:
        for _, group in duplicates.groupby(COLS["same_run"]):
            print("Duplicate rows found:")
            first_row = True
            for index, row in group.iterrows():
                if first_row:
                    print(
                        f"Instance: {row['instance_name']}, Policies: {row['policies']}, Budget: {row['budget']}, Solver: {row['solver']}"
                    )
                    first_row = False
                print(
                    "    --> "
                    + ", ".join(
                        [
                            f"{col}: {row[col]}"
                            for col in COLS["solution_outputs"] + ["time"]
                        ]
                    )
                )
            combined_df = combined_df.drop_duplicates(subset=match, keep="first")
    return combined_df


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate new raw run dataframe to main experiment dataframe."
    )
    parser.add_argument(
        "existing", metavar="X", type=str, nargs=1, help="existing dataframe to add to"
    )
    parser.add_argument(
        "new", metavar="N", type=str, nargs="+", help="new dataframe to add"
    )

    args = parser.parse_args()
    results_path = args.existing[0]
    results_df = pd.read_csv(results_path)

    add_cols = {
        key: -1 for key in COLS["processed"] if key not in set(results_df.columns)
    }
    for col, val in add_cols.items():
        results_df[col] = val

    for new_path in args.new:
        new_df = pd.read_csv(new_path)
        add_cols = {key: -1 for key in COLS["processed"] if key not in new_df.columns}
        for col, val in add_cols.items():
            new_df[col] = val
        results_df = concatenate_dataframes(results_df, new_df)

    results_df.to_csv(results_path, columns=COLS["processed"], header=True, index=False)


if __name__ == "__main__":
    main()
