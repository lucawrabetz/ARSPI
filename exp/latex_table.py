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

class OutputStructure:
    """
    Object to abstract misc params + columns about the output.
    Args:
        exp_type: experiment type (e.g. 'all_algorithms')
    Attributes:
        BE_ratio (bool)
        approximation_ratio (bool)
        sp_column (str)
        solvers (list[str])
        exact_solvers (list[str])
        time_columns (list[str])
        input_columns (list[str])
        int_columns (list[str])
        increment_columns (list[str])
        rational_columns (list[str])
        output_columns (list[str])
    """
    def __init__(self, exp_type):
        self.BE_ratio = False
        self.approximation_ratio = True
        if exp_type == MIPSYM_TYPENAME:
            self.approximation_ratio = False
        self.sp_column = 'MIP_objective'
        if exp_type == BE_TYPENAME:
            self.BE_ratio = True
            self.sp_column = 'BENDERS_objective'
        self.solvers = []
        if exp_type == ALL_TYPENAME:
            self.solvers = ['MIP', 'BENDERS', 'ENUMERATION', 'GREEDY']
        elif exp_type == MIPSYM_TYPENAME:
            self.solvers = ['MIP']
        elif exp_type == BE_TYPENAME:
            self.solvers = ['BENDERS', 'ENUMERATION', 'GREEDY']
        self.exact_solvers = [sol for sol in self.solvers if sol != 'GREEDY']
        self.time_columns = [sol + '_time' for sol in self.solvers]
        self.input_columns = ['k_zero', 'nodes', 'arcs', 'scenarios', 'policies']
        self.objective_columns = [sol + '_objective' for sol in self.solvers]
        self.int_columns = self.input_columns + self.objective_columns

        self.increment_columns = [sol + '_adaptiveincrement' for sol in self.solvers]
        self.rational_columns = self.increment_columns + self.time_columns
        MIP_OUT = ['MIP_time', 'MIP_gap', 'MIP_adaptiveincrement']
        BENDERS_OUT = ['BENDERS_time', 'BENDERS_gap',
                       'BENDERS_cuts_rounds', 'BENDERS_adaptiveincrement']
        ENUMERATION_OUT = ['ENUMERATION_time', 'ENUMERATION_adaptiveincrement']
        GREEDY_OUT = ['GREEDY_time', 'approximation_ratio', 'GREEDY_adaptiveincrement']
        self.output_columns = self.input_columns

        if 'MIP' in self.solvers:
            self.output_columns.extend(MIP_OUT)
            self.rational_columns.append('MIP_gap')
        if 'BENDERS' in self.solvers:
            self.output_columns.extend(BENDERS_OUT)
            self.rational_columns.append('BENDERS_cuts_rounds')
            self.rational_columns.append('BENDERS_gap')
        if 'ENUMERATION' in self.solvers:
            self.output_columns.extend(ENUMERATION_OUT)
        if self.BE_ratio:
            self.output_columns.append('BE_ratio')
            self.rational_columns.append('BE_ratio')
        if 'GREEDY' in self.solvers:
            self.rational_columns.append('approximation_ratio')
            self.output_columns.extend(GREEDY_OUT)


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


def write_latex_table(args, df, out):
    latex_csv = append_date('final') + '.csv'
    latex_csv_path = os.path.join('latex', latex_csv)
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
    parser = argparse.ArgumentParser(
        prog='AspiLatexTable',
        description='Create a latex table from a results csv for an Aspi computational experiment.')
    parser.add_argument('experiment_type', metavar='E', type=str, nargs=1, default=ALL_TYPENAME,
                        choices=[ALL_TYPENAME, BE_TYPENAME, MIPSYM_TYPENAME],
                        help='name of experiment type - ' + ALL_TYPENAME + ' or ' + BE_TYPENAME + ' or ' + MIPSYM_TYPENAME)
    args = parser.parse_args()
    run_df = read_results_data()
    out_structure = OutputStructure(args.experiment_type[0])
    # if out_structure.BE_ratio: compute_bendersenum_ratio(run_df)
    avg_df = group_by_instance_average_by_run(run_df)
    post_cleanup(avg_df, out_structure)
    latex_table_path = write_latex_table(args, avg_df, out_structure)
    final_table_cleanup(latex_table_path)


if __name__ == "__main__":
    main()
