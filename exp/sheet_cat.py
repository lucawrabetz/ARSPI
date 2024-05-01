from numpy import append
import pandas as pd
import shutil
import argparse

from lib.util import *
from lib.paths import FINALCSVPATH, BACKUPS

def backup_finalcsv():
    filename = check_make_dir(append_date(FINALCSVFILE.split(".")[0]), 0, False)
    shutil.copyfile(FINALCSVPATH, os.path.join(BACKUPS, filename))

def check_columns(A, B):
    if set(A.columns) != set(B.columns):
        print(A.columns)
        print(B.columns)
        raise ValueError(
            "Columns of DataFrame A do not match columns of DataFrame B. Cannot concatenate."
        )


def duplicate_check(df, cols):
    duplicates = df[df.duplicated(subset=cols, keep=False)]
    print("Full duplicates found: " + str(len(duplicates)) + "\n")
    if not duplicates.empty:
        print(
            "Found " + str(len(duplicates)) + " strongly matching rows, dropping.\n"
        )
        df.drop_duplicates(subset=cols, keep="first", inplace=True)


def concatenate_dataframes(A, B):
    check_columns(A, B)
    combined_df = pd.concat([A, B], ignore_index=True, sort=False)
    cols_names_outputs = [c.name for c in COLS["same_run_and_outputs"]]
    duplicate_check(combined_df, cols_names_outputs)
    return combined_df


def main():
    # The entities to be defined:
    # UI/Parser
    # DataFrameReader
    # DataFrameCleaner
    # DataFrameAppender
    # DataFrameWriter
    # CustomPlotSession (exists)

    # UI/Parser takes arguments from user and holds on to them.
    parser = argparse.ArgumentParser(
        description="Concatenate new raw run dataframe to main experiment dataframe."
    )
    parser.add_argument(
        "new", metavar="N", type=str, nargs="+", help="new dataframe to add"
    )
    args = parser.parse_args()
    backup_finalcsv()

    # Parser hands existing file path to Reader, reads and stores existing dataframe
    results_df = pd.read_csv(FINALCSVPATH)

    # Reader hands raw data to Cleaner (??)
    cleanup_to_processed(results_df)

    # Cleaner hands clean existing data to appender
    for new_path in args.new:
        # Parser hands new file path to Reader, reads and stores new dataframe
        new_df = pd.read_csv(new_path)
        # Reader hands raw data to ()
        cleanup_to_processed(new_df)
        # () hands clean new data to Appender, appends.
        results_df = concatenate_dataframes(results_df, new_df)

    # Appender hands final data to Writer, who writes.
    final_write(results_df, FINALCSVPATH)


if __name__ == "__main__":
    main()
