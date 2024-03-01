# Function + required dict to more generally map args to solvers.
# Case insensitive.
RAW_COLUMNS = ["set_name", "instance_name", "nodes", "arcs", "k_zero", "density", "scenarios", "budget", "policies", "solver", "unbounded", "optimal", "objective", "gap", "time", "cuts_rounds", "cuts_added", "avg_cbtime", "avg_sptime", "partition", "m_sym", "g_sym"]

RATIO_COLUMNS = ["best_objective", "best_optimal", "empirical_suboptimal_ratio", "empirical_optimal_ratio"]

SAME_RUN_COLUMNS = ["set_name", "instance_name", "nodes", "arcs", "k_zero", "density", "scenarios", "budget", "policies", "solver"]

TIME_COLUMNS = ["time", "avg_cbtime", "avg_sptime"]
SECONDS_COLUMNS = [i + "_s" for i in TIME_COLUMNS]

INT_COLUMNS = ["nodes", "arcs", "k_zero", "scenarios", "budget", "policies", "cuts_rounds", "cuts_added", "m_sym", "g_sym"]

RATIONAL_COLUMNS = ["objective", "gap", "density"] + TIME_COLUMNS + SECONDS_COLUMNS + RATIO_COLUMNS

DP = {
    "objective": 2,
    "gap": 2,
    "density": 2,
    "time": 2,
    "avg_cbtime": 2,
    "avg_sptime": 2,
    "time_s": 2,
    "avg_cbtime_s": 2,
    "avg_sptime_s": 2,
    'best_objective': 2,
    'best_optimal': 2,
    'empirical_suboptimal_ratio': 2,
    'empirical_optimal_ratio': 2
}

FINAL_COLUMNS = RAW_COLUMNS + SECONDS_COLUMNS + RATIO_COLUMNS

SOLVER_FLAGS = {
    "MIP" : set({"m", "mip", "sp"}),
    "BENDERS" : set({"b", "benders",}),
    "ENUMERATION" : set({"e", "enum", "enumeration"}),
    "GREEDY" : set({"g", "greedy"}),
}

FINAL_COLUMN_TO_ITSELF = {
     'set_name': 'set_name',
    'instance_name': 'instance_name',
    'nodes': 'nodes',
    'arcs': 'arcs',
    'k_zero': 'k_zero',
    'density': 'density',
    'scenarios': 'scenarios',
    'budget': 'budget',
    'policies': 'policies',
    'solver': 'solver',
    'unbounded': 'unbounded',
    'optimal': 'optimal',
    'objective': 'objective',
    'gap': 'gap',
    'time': 'time',
    'cuts_rounds': 'cuts_rounds',
    'cuts_added': 'cuts_added',
    'avg_cbtime': 'avg_cbtime',
    'avg_sptime': 'avg_sptime',
    'partition': 'partition',
    'm_sym': 'm_sym',
    'g_sym': 'g_sym',
    'time_s': 'time_s', 
    'avg_cbtime_s': 'time_s', 
    'avg_sptime_s': 'time_s', 
    'best_objective': 'best_objective',
    'best_optimal': 'best_optimal',
    'empirical_suboptimal_ratio': 'empirical_suboptimal_ratio',
    'empirical_optimal_ratio': 'empirical_optimal_ratio'
}

FINAL_COLUMN_TO_PRETTY_PRINT = {
    'set_name': 'Set Name', 
    'instance_name': 'Instance Name', 
    'nodes': 'Nodes', 
    'arcs': 'Arcs', 
    'k_zero': 'Groups', 
    'density': 'Density', 
    'scenarios': 'Followers', 
    'budget': 'Budget', 
    'policies': 'Policies', 
    'solver': 'Solver', 
    'unbounded': 'Unbounded', 
    'optimal': 'Optimal', 
    'objective': 'Objective', 
    'gap': 'Gap', 
    'time': 'Running Time (ms)', 
    'cuts_rounds': 'Callbacks', 
    'cuts_added': 'Cuts Added', 
    'avg_cbtime': 'Average Callback Time (ms)', 
    'avg_sptime': 'Average Cut Separation Time (ms)', 
    'partition': 'Partition', 
    'm_sym': 'Manual Symmetry', 
    'g_sym': 'Gurobi Symmetry',
    'time_s': 'Running Time (s)', 
    'avg_cbtime_s': 'Average Callback Time (s)', 
    'avg_sptime_s': 'Average Cut Separation Time (s)', 
    'best_objective': 'Best Objective',
    'best_optimal': 'Best Optimal Objective',
    'empirical_suboptimal_ratio': 'Empirical Suboptimal Ratio',
    'empirical_optimal_ratio': 'Empirical Optimal Ratio'
}

FINAL_COLUMN_TO_SHORT_PRINT = {
    'set_name': 'Set', 
    'instance_name': 'Instance', 
    'nodes': 'N', 
    'arcs': 'M', 
    'k_zero': 'k0', 
    'density': 'D', 
    'scenarios': 'P', 
    'budget': 'r0', 
    'policies': 'K', 
    'solver': 'Sol', 
    'unbounded': 'Unb', 
    'optimal': 'Opt', 
    'objective': 'Obj', 
    'gap': 'Gap (%)', 
    'time': 'T (ms)', 
    'cuts_rounds': 'Cb', 
    'cuts_added': 'Cuts', 
    'avg_cbtime': 'CbT (ms)', 
    'avg_sptime': 'CutT (ms)', 
    'partition': 'Part', 
    'm_sym': 'Msym', 
    'g_sym': 'Gsym',
    'time_s': 'T (s)', 
    'avg_cbtime_s': 'CbT (s)', 
    'avg_sptime_s': 'CutT (s)', 
    'best_objective': 'MaxObj',
    'best_optimal': 'MaxOpt',
    'empirical_suboptimal_ratio': 'rObj',
    'empirical_optimal_ratio': 'rOpt'
}

def get_solver_from_flag(flag):
    for solver, flag_set in SOLVER_FLAGS.items():
        if flag.lower() in flag_set: return solver
    return None



def ms_to_s(df, col):
    '''
    Convert column col in df from ms to s.
    '''
    new_col = col + "_s"
    df[new_col] = df[col].div(1000)

def add_seconds_columns(df):
    '''
    Convert all time columns to s.
    '''
    for col in TIME_COLUMNS:
        ms_to_s(df, col)

def round_int_columns(df):
    '''
    Round all integer columns.
    '''
    for col in INT_COLUMNS:
        df[col] = df[col].astype(int)

def round_dp_columns(df):
    '''
    Round all decimal columns to DP decimal points.
    '''
    for col in RATIONAL_COLUMNS:
        df[col] = df[col].astype(float)
        df[col] = df[col].round(DP[col])

def round(df):
    round_int_columns(df)
    round_dp_columns(df)

def main():
    pass

if __name__=='__main__':
    for i in RAW_COLUMNS:
        print("'" + i + "': " + "'" + i + "'")

    main()
