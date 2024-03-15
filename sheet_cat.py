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


def print_duplicate_group(group):
    print("Duplicate rows found to inspect:")
    first_row = True
    for _, row in group.iterrows():
        if first_row:
            print_header(row)
            first_row = False
        print_finished_row(row)


def duplicate_check(df, cols, drop=False):
    # Will only modify df if drop=True, otherwise the user will
    # get an option to exit the program.
    duplicates = df[df.duplicated(subset=cols, keep=False)]
    if not duplicates.empty:
        if drop:
            print(
                "Found " + str(len(duplicates)) + " strongly matching rows, dropping.\n"
            )
            df.drop_duplicates(subset=cols, keep="first", inplace=True)
        else:
            for _, group in duplicates.groupby(cols):
                print_duplicate_group(group)
                # TODO: consider expanding continue into keep or drop.
                user_input = input("Press 'c' to continue or 'e' to exit: ")
                if user_input.lower() == "e":
                    raise KeyboardInterrupt("Exiting due to duplicates.")


def concatenate_dataframes(A, B):
    check_columns(A, B)
    combined_df = pd.concat([A, B], ignore_index=True, sort=False)
    cols_names_outputs = [c.name for c in COLS["same_run_and_outputs"]]
    cols_names_params = [c.name for c in COLS["same_run_params"]]
    duplicate_check(combined_df, cols_names_outputs, drop=True)
    duplicate_check(combined_df, cols_names_params)
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

    common_cleanup(results_df)
    for new_path in args.new:
        new_df = pd.read_csv(new_path)
        common_cleanup(new_df)
        results_df = concatenate_dataframes(results_df, new_df)

    final_write(results_df, results_path)


if __name__ == "__main__":
    main()
