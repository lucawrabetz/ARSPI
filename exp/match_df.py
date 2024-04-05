import pandas as pd
import argparse

def read_dataframe(filepath):
    try:
        df = pd.read_csv(filepath)
        return df
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return None

def read_dataframes(filepath1, filepath2):
    df1 = read_dataframe(filepath1)
    df2 = read_dataframe(filepath2)

    if df1 is None or df2 is None:
        print("Error: One or both dataframes could not be loaded.")
        return None, None

    if df1.shape != df2.shape:
        raise ValueError("Dataframes have different dimensions. Cannot compare.")

    return df1, df2

def compare_dataframes(df1, df2):
    if df1.equals(df2):
        return "MATCH"
    else:
        mismatch_info = []
        diff_locations = [(i, j) for i in range(df1.shape[0]) for j in range(df1.shape[1]) if df1.iloc[i, j] != df2.iloc[i, j]]
        if diff_locations:
            diff_info = [f"[{i}, {df1.columns[j]}]: {df1.iloc[i, j]} != {df2.iloc[i, j]}" for i, j in diff_locations]
            mismatch_info.append("Differences found at:\n" + "\n".join(diff_info))
        return "\n".join(mismatch_info)

def main():
    parser = argparse.ArgumentParser(description="Compare two dataframes from filepaths")
    parser.add_argument("file1", type=str, help="Path to the first file")
    parser.add_argument("file2", type=str, help="Path to the second file")
    args = parser.parse_args()
    df1, df2 = read_dataframes(args.file1, args.file2)
    result = compare_dataframes(df1, df2)
    print(result)

if __name__ == "__main__":
    main()
