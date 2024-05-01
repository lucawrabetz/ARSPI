import os
import argparse
import pandas as pd
from lib.util import *

# Names of experiment instance sets.
ALL_INSTANCES = 'experiment_allalgorithms'
MIPSYM_INSTANCES = 'experiment_allalgorithms_symmetrynondecreasinggurobiaggressive'
BE_INSTANCES = 'experiment_bendersenum'
LAYERRATIO_INSTANCES = 'experiment_layerratio'

# Names of experiment output types.
ALL_TYPENAME = 'all_algorithms'
BE_TYPENAME = 'benders_enum'
MIPSYM_TYPENAME = 'mip_symmetry'

# Decimal points to be rounded to for decimal columns.
DP = 2

def compute_bendersenum_ratio(df):
    '''
    Compute new column 'BE_ratio' reflecting the ratio in objective values Benders / Enumeration.
    '''
    df['BE_ratio'] = df['BENDERS_objective'] / df['ENUMERATION_objective']


def average_by_run(group):
    new_rows = []
    run_groups = group.groupby(['instance_name', 'budget', 'policies', 'subsolver', 'm_sym', 'g_sym'])
    for name, group in run_groups:
        policies = name[2]
        if policies == 0:
            continue
        mean_df = group.mean(numeric_only=True)
        new_rows.append(mean_df)
    return new_rows


def group_by_instance_average_by_run(df):
    '''
    Group data frame by instance type ['set_name', 'nodes', 'k_zero', 'scenarios'].
                                        - also implies that ['arcs', 'density'] will match. (INSTANCE PARAMETERS -> INPUT -> FEATURE).
    For each group, the data can differ in the input columns : ['instance_name', 'budget', 'policies', 'subsolver', 'm_sym', 'g_sym']. (RUN PARAMETERS -> INPUT -> FEATURE).
    For each group, the outputs can all differ. (OUTPUT -> FEATURE).
    Average the outputs (10 runs) for each run group (average_by_run(group)).
    '''
    instance_groups = df.groupby(['set_name', 'nodes', 'k_zero', 'scenarios'])
    avg_rows = []
    for name, group in instance_groups:
        # name is a tuple of the form (set_name, nodes, k_zero, scenarios)
        # max_k = grouped['policies'].max()
        new_rows = average_by_run(group)
        avg_rows.extend(new_rows)
    avg_df = pd.DataFrame(avg_rows)
    return avg_df


def read_results_data():
    csv = "final.csv"
    df = pd.read_csv(csv)
    return df


def write_latex_table(df):
    latex_csv = append_date('final') + '.csv'
    latex_csv_path = os.path.join('latex', latex_csv)
    # TODO: replace out.output_columns with every solver's 
    # output columns.
    df.to_csv(path_or_buf=latex_csv_path, sep='&',
              columns=out.output_columns, index=False)
    return latex_csv_path


def final_table_cleanup(path):
    output_path = path.split('.')[0] + '-clean.csv'

    with open(path, 'r') as input_file, open(output_path, 'w') as output_file:
        # Read and write the header without modification
        header = input_file.readline().strip()
        output_file.write(header + '\n')
        # Process the rest of the lines
        for line in input_file:
            # Strip any leading or trailing whitespace, and append '\\'
            modified_line = line.strip() + '\\\\'
            output_file.write(modified_line + '\n')

def main():
    parser = FeatureArgParser("Generate Latex Source for table from filtered data.")
    parser.add_feature_args(COLS["processed"])
    run_df = read_results_data()
    # filtering
    avg_df = group_by_instance_average_by_run(run_df)
    latex_table_path = write_latex_table(args, avg_df, out_structure)
    final_table_cleanup(latex_table_path)


if __name__ == "__main__":
    main()
