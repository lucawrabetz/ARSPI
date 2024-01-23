import pandas as pd
import argparse


def concatenate_dataframes(A, B):
    if not A.columns.equals(B.columns):
        raise ValueError(
            "Columns of DataFrame A do not exactly match columns of DataFrame B. Cannot concatenate.")
    result = pd.concat([A, B], ignore_index=True)
    return result


def main():
    parser = argparse.ArgumentParser(
        description='Concatenate new raw run dataframe to main experiment dataframe.')
    parser.add_argument('existing', metavar='X', type=str, nargs=1,
                        help='existing dataframe to add to')
    parser.add_argument('new', metavar='N', type=str, nargs='+',
                        help='new dataframe to add')

    args = parser.parse_args()
    results_path = args.existing[0]
    results_df = pd.read_csv(results_path)

    for new_path in args.new:
        new_df = pd.read_csv(new_path)
        results_df = concatenate_dataframes(results_df, new_df)

    results_df.to_csv(results_path, header=True, index=False)


if __name__ == '__main__':
    main()
