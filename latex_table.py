import os
import argparse
import pandas as pd

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


def ms_to_s(df, col):
    '''
    Convert column col in df from ms to s.
    '''
    df[col] = df[col].div(1000)


def round_int_columns(df, out):
    '''
    Round all integer columns.
    '''
    for col in out.int_columns:
        df[col] = df[col].astype(int)


def convert_all_times(df, out):
    '''
    Convert all time columns to s.
    '''
    for col in out.time_columns:
        ms_to_s(df, col)


def round_dp_columns(df, out):
    '''
    Round all decimal columns to DP decimal points.
    '''
    for col in out.rational_columns:
        df[col] = df[col].astype(float)
        df[col] = df[col].round(DP)


def post_cleanup(df, out):
    '''
    Complete any final rounding, units conversions.
    '''
    round_int_columns(df, out)
    convert_all_times(df, out)
    round_dp_columns(df, out)


def compute_approximation_ratio(df, out):
    '''
    Compute empirical approximation ratio using an exact solution.
    '''
    exact_columns = [solver + '_objective' for solver in out.exact_solvers]
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


def compute_adaptive_increment(df, out):
    '''
    Compute the 'adaptive increment', that is, the ratio between the objective value of a solution and the uninterdicted shortest path for that instance.
    Compute this for every solver's output as passed in the parameter solvers, and add it as a new column solver+'_adaptiveincrement' for each.
    '''
    add_shortest_path_column(df, out.sp_column)
    for solver in out.solvers:
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

def write_latex_table(args, df, out):
    latex_csv = args.set_name[0] + '-latex.csv'
    latex_csv_path = os.path.join('results', latex_csv)
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
    parser.add_argument('set_name', metavar='S', type=str, nargs=1, default=ALL_INSTANCES,
                        choices=[ALL_INSTANCES, BE_INSTANCES, MIPSYM_INSTANCES, LAYERRATIO_INSTANCES],
                        help='name of experiment set - ' + ALL_INSTANCES + ' or ' + BE_INSTANCES + ' or ' + LAYERRATIO_INSTANCES + ' or '+ MIPSYM_INSTANCES)
    parser.add_argument('experiment_type', metavar='E', type=str, nargs=1, default=ALL_TYPENAME,
                        choices=[ALL_TYPENAME, BE_TYPENAME, MIPSYM_TYPENAME],
                        help='name of experiment type - ' + ALL_TYPENAME + ' or ' + BE_TYPENAME + ' or ' + MIPSYM_TYPENAME)
    args = parser.parse_args()
    run_df = read_results_data(args)
    out_structure = OutputStructure(args.experiment_type[0])
    if out_structure.approximation_ratio: compute_approximation_ratio(run_df, out_structure)
    compute_adaptive_increment(run_df, out_structure)
    if out_structure.BE_ratio: compute_bendersenum_ratio(run_df)
    avg_df = take_averages(run_df)
    post_cleanup(avg_df, out_structure)
    latex_table_path = write_latex_table(args, avg_df, out_structure)
    final_table_cleanup(latex_table_path)


if __name__ == "__main__":
    main()
