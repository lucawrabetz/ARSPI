import os
import argparse
import pandas as pd

# Decimal points to be rounded to for decimal columns.
DP = 2
BE_ratio = True
EXACT_SOLVERS = [
    #    'MIP',
    'BENDERS',
    'ENUMERATION'
]
APPROXIMATE_SOLVERS = [
    #    'GREEDY'
]
SOLVERS = EXACT_SOLVERS + APPROXIMATE_SOLVERS

TIME_COLUMNS = [solver + '_time' for solver in SOLVERS]
PARAM_COLUMNS = ['k_zero', 'nodes', 'arcs', 'scenarios', 'policies']
OBJECTIVE_COLUMNS = [solver + '_objective' for solver in SOLVERS]
INT_COLUMNS = PARAM_COLUMNS + OBJECTIVE_COLUMNS

INCREMENT_COLUMNS = [solver + '_adaptiveincrement' for solver in SOLVERS]
DECIMAL_COLUMNS = INCREMENT_COLUMNS + TIME_COLUMNS

MIP_OUT = ['MIP_time', 'MIP_gap', 'MIP_adaptiveincrement']
BENDERS_OUT = ['BENDERS_time', 'BENDERS_gap',
               'BENDERS_cuts_rounds', 'BENDERS_adaptiveincrement']
ENUMERATION_OUT = ['ENUMERATION_time', 'ENUMERATION_adaptiveincrement']
GREEDY_OUT = ['GREEDY_time', 'approximation_ratio', 'GREEDY_adaptiveincrement']
OUTPUT_COLUMNS = PARAM_COLUMNS

if 'MIP' in SOLVERS:
    OUTPUT_COLUMNS.extend(MIP_OUT)
    DECIMAL_COLUMNS.append('MIP_gap')
if 'BENDERS' in SOLVERS:
    OUTPUT_COLUMNS.extend(BENDERS_OUT)
    DECIMAL_COLUMNS.append('BENDERS_cuts_rounds')
    DECIMAL_COLUMNS.append('BENDERS_gap')
if 'ENUMERATION' in SOLVERS:
    OUTPUT_COLUMNS.extend(ENUMERATION_OUT)
if 'GREEDY' in SOLVERS:
    DECIMAL_COLUMNS.append('approximation_ratio')
    OUTPUT_COLUMNS.extend(GREEDY_OUT)
if BE_ratio:
    DECIMAL_COLUMNS.append('BE_ratio')
    OUTPUT_COLUMNS.append('BE_ratio')


def ms_to_s(df, col):
    '''
    Convert column col in df from ms to s.
    '''
    df[col] = df[col].div(1000)


def round_int_columns(df):
    '''
    Round all integer columns.
    '''
    for col in INT_COLUMNS:
        df[col] = df[col].astype(int)


def convert_all_times(df):
    '''
    Convert all time columns to s.
    '''
    for col in TIME_COLUMNS:
        ms_to_s(df, col)


def round_dp_columns(df):
    '''
    Round all decimal columns to DP decimal points.
    '''
    for col in DECIMAL_COLUMNS:
        df[col] = df[col].astype(float)
        df[col] = df[col].round(DP)


def post_cleanup(df):
    '''
    Complete any final rounding, units conversions.
    '''
    round_int_columns(df)
    convert_all_times(df)
    round_dp_columns(df)


def compute_approximation_ratio(df):
    '''
    Compute empirical approximation ratio using an exact solution.
    '''
    exact_columns = [solver + '_objective' for solver in EXACT_SOLVERS]
    max_exact_objective = df[exact_columns].max(axis=1)
    df['approximation_ratio'] = df['GREEDY_objective'] / max_exact_objective


def add_shortest_path_column(df, objective_col):
    '''
    Add a temporary column with the uninterdicted shortest path for every instance.
    '''
    # Find the MIP_objective value where policies = 0 for each instance_name (returning a DataFrame with columns 'instance_name' and 'MIP_objective').
    shortest_paths = df.loc[df['policies'] ==
                            0, ['instance_name', objective_col]]
    # Create a dictionary mapping instance_name to its corresponding uninterdicted shortest path.
    shortest_paths_dict = shortest_paths.set_index(
        'instance_name')[objective_col].to_dict()
    df['shortest_path'] = df['instance_name'].map(shortest_paths_dict)


def compute_adaptive_increment(df, objective_col='MIP_objective'):
    '''
    Compute the 'adaptive increment', that is, the ratio between the objective value of a solution and the uninterdicted shortest path for that instance.
    Compute this for every solver's output as passed in the parameter solvers, and add it as a new column solver+'_adaptiveincrement' for each.
    '''
    add_shortest_path_column(df, objective_col)
    for solver in SOLVERS:
        newcolumnname = solver + '_adaptiveincrement'
        objectivecolumnname = solver + '_objective'
        df[newcolumnname] = df[objectivecolumnname] / df['shortest_path']


def compute_bendersenum_ratio(df):
    '''
    Compute new column 'BE_ratio' reflecting the ratio in objective values Benders / Enumeration.
    '''
    df['BE_ratio'] = df['BENDERS_objective'] / df['ENUMERATION_objective']


def k_meanrows(group):
    '''
    Return K rows averaging the 10 runs for each value of k in the DataFrame 'group'.
    '''
    new_rows = []
    grouped_by_k = group.groupby('policies')
    for k, df in grouped_by_k:
        if k == 0:
            continue
        mean_df = df.mean(numeric_only=True)
        new_rows.append(mean_df)
    return new_rows


def take_averages(df):
    '''
    For every instance 'type' (matching nodes, k_zero, scenarios), average the 10 instance ids into K new rows,
    where K is the number of different policies ran for each instance in the instance type set.
    '''
    grouped = df.groupby(['nodes', 'k_zero', 'scenarios'])
    avg_rows = []
    for name, group in grouped:
        # max_k = grouped['policies'].max()
        new_rows = k_meanrows(group)
        avg_rows.extend(new_rows)
    avg_df = pd.DataFrame(avg_rows)
    return avg_df


def read_results_data(args):
    csv = args.set_name[0] + '.csv'
    csvpath = os.path.join('results', csv)
    df = pd.read_csv(csvpath)
    return df


def write_latex_table(args, df):
    latex_csv = args.set_name[0] + '-latex.csv'
    latex_csv_path = os.path.join('results', latex_csv)
    df.to_csv(path_or_buf=latex_csv_path, sep='&',
              columns=OUTPUT_COLUMNS, index=False)


# ALLCOLUMNS = ["set_name", "instance_name", "nodes", "arcs", "k_zero", "density", "scenarios", "budget", "policies",
#               "MIP_OPTIMAL", "MIP_objective", "MIP_gap", "MIP_time", "MIP_partition",
#               "BENDERS_OPTIMAL", "BENDERS_objective", "BENDERS_gap", "BENDERS_time", "BENDERS_cuts_rounds", "BENDERS_partition",
#               "ENUMERATION_OPTIMAL", "ENUMERATION_objective", "ENUMERATION_time", "ENUMERATION_partition",
#               "GREEDY_objective", "GREEDY_time", "GREEDY_partition"
#               ]
#
# # output columns
# COLUMNS = ["k_zero",
#            "nodes", "arcs",
#            "scenarios", "policies",
#            "MIP_objective", "MIP_time", "MIP_gap",
#            "BENDERS_objective", "BENDERS_time", "BENDERS_gap", "BENDERS_cuts_rounds",
#            "ENUMERATION_objective", "ENUMERATION_time",
#            "GREEDY_objective", "GREEDY_time", "approximation_ratio"]
#
# NUM_COLUMNS = ["MIP_objective", "MIP_time", "MIP_gap", "MIP_adaptiveincrement",
#                "BENDERS_objective", "BENDERS_time", "BENDERS_gap", "BENDERS_cuts_rounds", "BENDERS_adaptiveincrement",
#                "ENUMERATION_objective", "ENUMERATION_time", "ENUMERATION_adaptiveincrement",
#                "GREEDY_objective", "GREEDY_time", "approximation_ratio", "GREEDY_adaptiveincrement"]
#
# # output columns for averages
# AVG_COLUMNS = ["k_zero",
#                "nodes", "arcs",
#                "scenarios", "policies",
#                "MIP_time", "MIP_gap", "MIP_adaptiveincrement",
#                "BENDERS_time", "BENDERS_gap", "BENDERS_cuts_rounds", "BENDERS_adaptiveincrement",
#                "ENUMERATION_time", "ENUMERATION_adaptiveincrement",
#                "GREEDY_time", "approximation_ratio", "GREEDY_adaptiveincrement"]


# def convert_numerical_to_float(data):
#     for col in NUM_COLUMNS:
#         data[col] = data[col].astype(float)


def main():
    parser = argparse.ArgumentParser(
        description='Create a latex table from a results csv.')
    parser.add_argument('set_name', metavar='S', type=str, nargs=1,
                        help='name of experiment')
    args = parser.parse_args()
    run_df = read_results_data(args)
    # run_df.replace(",", "", regex=True, inplace=True)
    # compute_approximation_ratio(run_df)
    compute_adaptive_increment(run_df, objective_col='BENDERS_objective')
    compute_bendersenum_ratio(run_df)
    # convert_numerical_to_float(run_df)
    avg_df = take_averages(run_df)
    post_cleanup(avg_df)

    write_latex_table(args, avg_df)


if __name__ == "__main__":
    main()
