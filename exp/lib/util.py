import os
import pandas as pd
import logging
from typing import Dict, Any
from lib.feature import COLS

DP = 2

def check_make_dir(path: str, i: int, makedir: bool = True):
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
        if makedir:
            os.mkdir(path + "-" + str(i))
        return path + "-" + str(i)

def ms_to_s(df: pd.DataFrame, col: str):
    """
    Convert column col in df from ms to s.
    """
    new_col = col + "_s"
    if df[col] is None:
        logging.warning(f"Input time ms col {col} is None, conversion to seconds aborted.")
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
            df[col] = df[col].round(DP)

def round(df: pd.DataFrame):
    round_int_columns(df)
    round_dp_columns(df)

def print_dict(d: Dict[Any, Any]):
    print("{")
    for k, v in d.items():
        print(k + ": " + v + ",")
    print("}")

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

def remove_samerun_duplicates(df):
    """
    Remove duplicates for the same exact run parameters and hyperparameters.
    Keep the row with the highest objective, splitting ties by running time (minimum chosen).
    Make changes to df in place.
    """
    matching_colnames = [f.name for f in COLS["same_run_hyperparams"]]
    df.sort_values(
        by=["objective", "time"], ascending=[False, True], inplace=True
    )
    df.drop_duplicates(subset=matching_colnames, keep="first", inplace=True)
    logging.info("Removed duplicates for the same run.")

def cleanup_to_processed(df):
    """
    Add any missing columns using defaults to meet COLS["processed"].
    """
    add_cols = {
        key.name: key.default
        for key in COLS["processed"]
        if key.name not in set(df.columns)
    }
    for col, val in add_cols.items():
        print("Adding column ", col, " with default ", val, ".")
        df[col] = val
    pass

def cleanup_to_finished(df) -> Any:
    """
    Remove same run dups, and add quick processing columns to meet COLS["finished"].
    Returns a new dataframe to use, so as not to make changes to raw df.
    """
    data_df = df.copy()
    remove_samerun_duplicates(df)
    add_seconds_columns(data_df)
    return data_df

def cleanup_to_pretty(df) -> Any:
    """
    Rounding and stuff that only matters for printing / tables etc, but we don't want to do the raw data.
    """
    pretty_df = df.copy()
    round(pretty_df)
    return pretty_df

def final_write(df, path):
    """
    Always starts with a df with "finished" columns, because common_clenup was called.
    """
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
