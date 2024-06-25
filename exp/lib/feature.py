import argparse
import pandas as pd
import logging
from typing import Dict, Any, List, Tuple, Set, Type, TypeVar, Callable

FEATURE_TYPE = TypeVar('FEATURE_TYPE')
# TODO: substitute the generic TypeVar FEATURE_TYPE with a list of allowed types for features.
# Consider doing this while also adding subclasses for some types of features.
# Check conversation "explicit type hinting:..." in my chat gpt channel.
class Feature:
    def __init__(self, name: str, feature_type: Type[FEATURE_TYPE], default: FEATURE_TYPE = None, pretty_output_name: str = None, compressed_output_name: str = None, allowed_values: List[FEATURE_TYPE] = []) -> None:
        self.name = name
        self.default = default
        self.type = feature_type
        if pretty_output_name: self.pretty_output_name = pretty_output_name
        else: self.pretty_output_name = name
        if compressed_output_name: self.compressed_output_name = compressed_output_name
        else: self.compressed_output_name = name
        if allowed_values: self.allowed_values = allowed_values
        else: self.allowed_values = []
        
        if type(default) != feature_type:
            logging.warning("Feature " + name + " - The type of default value " + str(default) + " does not match the specified feature type " + str(feature_type) + ".")

class Solver:
    def __init__(self, name: str):
        logging.info("Creating solver: " + name)
        self.name = name
        self.parameters: List[Feature] = []
        self.latex_output_features: List[Feature] = []
        self.commandline_flags: Set[str] = set()

    def add_parameters(self, parameters: List[Feature]):
        self.parameters.extend(parameters)

    def add_latex_output_features(self, latex_output_features: List[Feature]):
        self.latex_output_features.extend(latex_output_features)

    def add_commandline_flags(self, commandline_flags: Set[str]):
        self.commandline_flags.update(commandline_flags)

# TODO: should this function be a method of the Solver class? or of the FeatureArgParser class? Every way I spin it, it feels like a codesmell...
def flag_to_solver(flag: str) -> Solver:
    for solver in SOLVERS:
        if flag.lower() in solver.commandline_flags:
            return solver
    logging.warning("Solver not found for flag: " + flag + ".")
    return NONE

# TODO: split up this file into two files, one for the user to declare their solvers and features.
MIP = Solver("MIP")
MIP.add_commandline_flags({"m", "mip", "sp"})
BENDERS = Solver("BENDERS")
BENDERS.add_commandline_flags({"b", "benders"})
ENUMERATION = Solver("ENUMERATION")
ENUMERATION.add_commandline_flags({"e", "enum", "enumeration"})
GREEDY = Solver("GREEDY")
GREEDY.add_commandline_flags({"g", "greedy"})
NONE = Solver("NONE")
NONE.add_commandline_flags({"n", "none", "na"})

SOLVERS = [MIP, BENDERS, ENUMERATION, GREEDY, NONE]

class TypeMapper:
    """
    Class to map allowed types in the data model to string representations and vice versa.
    """
    def __init__(self):
        self.type_to_string: Dict[Type[Any], str] = {
            int: 'int',
            float: 'float',
            str: 'str',
            bool: 'bool',
            Solver: 'Solver',
        }
        self.string_to_type: Dict[str, Type[Any]] = {
            'int': int,
            'float': float,
            'str': str,
            'bool': bool,
            'Solver': Solver,
        }
        self.type_to_conversion = {
            int: int,
            float: float,
            str: lambda x: x,
            bool: self.convert_bool,
            Solver: self.convert_solver_type,
        }

    def convert_bool(self, default: str) -> bool:
        if default.lower() == 'true':
            return True
        elif default.lower() == 'false':
            return False
        else:
            raise ValueError(f"Default value {default} is not a valid boolean")

    def convert_solver_type(self, default: str) -> Solver:
        for s in SOLVERS:
            if default == s.name:
                return s
        logging.warning(f"No solver matches string value {default}.")
        return NONE

    def convert_default(self, default: str, data_type: Type[Any]) -> Any:
        """
        Convert the string default value to the same value of the correct matching type.
        """
        return self.type_to_conversion[data_type](default)

    def str_allowed(self, data_type: str) -> bool:
        """
        Check if a string is a valid data type.
        """
        return data_type in self.string_to_type.keys()

    def type_allowed(self, data_type: Type[Any]) -> bool:
        """
        Check if a type is a valid data type.
        """
        return data_type in self.type_to_string.keys()

    def type_to_str(self, data_type: Type[Any]) -> str:
        """
        Convert a type to its string representation.
        """
        return self.type_to_string[data_type]

    def str_to_type(self, data_type: str) -> Type[Any]:
        """
        Convert a string to its type representation.
        """
        return self.string_to_type[data_type]

ALLOWED_TYPES = TypeMapper()

class FeatureArgParser:
    '''
    Class to handle argument parsing where arguments are Feature - Value pairs.
    '''
    def __init__(self, description: str = "") -> None:
        self.parser = argparse.ArgumentParser(description=description)

    def construct_flag_strings(self, feature: Feature) -> Tuple[str, str, str]:
        arg_flag = "--" + feature.name
        help_msg = "Filter by " + feature.name + "."
        success_msg = "Added argument: " + arg_flag + " with type " + ALLOWED_TYPES.type_to_str(feature.type) + "."
        return (arg_flag, help_msg, success_msg)

    def add_solverfeature_arg(self, feature: Feature) -> str:
        (arg_flag, help_msg, success_msg) = self.construct_flag_strings(feature)
        self.parser.add_argument(arg_flag, type=flag_to_solver, choices=SOLVERS, help=help_msg)
        return success_msg

    def add_feature_arg(self, feature: Feature) -> str:
        (arg_flag, help_msg, success_msg) = self.construct_flag_strings(feature)
        self.parser.add_argument(arg_flag, type=feature.type, help=help_msg)
        return success_msg

    def add_storetruearg(self, arg_flag: str, help_msg: str = "") -> None:
        self.parser.add_argument(arg_flag, action="store_true", help=help_msg)

    def add_feature_args(self, features: List[Feature]) -> None:
        for f in features:
            if f.type == Solver:
                success_msg = self.add_solverfeature_arg(f)
            else:
                success_msg = self.add_feature_arg(f)
            logging.info(success_msg)

    def parse_args(self) -> argparse.Namespace:
        return self.parser.parse_args()

class DataFilterer:
    '''
    Class to handle filtering of main dataframe based on the Feature - value pairs held in a FeatureArgParser object.
    '''
    def __init__(self, args: argparse.Namespace) -> None:
        self.solver: Solver = args.solver
        self.args: Dict[Any, Any] = vars(args)

    def filter(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Filter dataframe based on provided kwargs.
        '''
        mask = pd.Series(True, index=df.index)
        # TODO: review this function now that we have a Solver class.
        for key, value in self.args.items():
            if (
                key != "file_path"
                and key != "average"
                and key != "solver"
                and key != "verbose"
                and value is not None
            ):
                mask = mask & (df[key] == value)
        if self.solver:
            mask = mask & (df["solver"] == self.solver.name)
        return df[mask]


# TODO: whether any of these "groups" need to be subclasses of features, or if they could be noted by some member flags. If that gets figured out, quite a few functions that take a list of features could iterate over all the features, and check their group belonging to decide whether to consider them, or how to treat them. Then we can remove all the lists and stuff.
# For the final interface:
# Each of the groups that follow have a single feature that is UNIQUE for every combination of the group.
# We call this feature the GROUP_ID, and it is always (IMPLEMENTATION) REQUIRED.
# There may be other required features in the group.
# We can always generate the GROUP_ID from the other features.
# As follows, group_name -> group_id
# INSTANCE FEATURES -> INSTANCE_NAME
# SOLVER -> SOLVER (This is a special group that should be a singleton) (TODO: not sure about this statement yet).
# SOLVER PARAMETERS -> PARAMETERS_ID
# TODO: for now only instance features have an ID (INSTANCE_NAME).
# TODO: all id_features should only be added as necessary, and they should be "finished" feature (they can be computed quickly), not stored in the database. This also means we should remove, e.g. INSTANCE_NAME from the "processed" columns and the database. Same goes for density, even though it is not an id_feature.
# TODO: id_generator -> instance_name(set_name, <all_other_inputs) -> instance_name.

# Instance Features - they define the INSTANCE OF THE PROBLEM.
# ID: INSTANCE_NAME
# REQUIRED: SET_NAME
SET_NAME = Feature("set_name", str, "layer", "Set Name", "Set")
INSTANCE_NAME = Feature("instance_name", str, "layer_123", "Instance Name", "Instance")
NODES = Feature("nodes", int, 123, "Nodes", "N")
ARCS = Feature("arcs", int, 1231, "Arcs", "M")
K_ZERO = Feature("k_zero", int, 5, "Groups", "k0")
DENSITY = Feature("density", float, 1231 / (123 * 122), "Density", "D")
SCENARIOS = Feature("scenarios", int, 5, "Followers", "P")
BUDGET = Feature("budget", int, 5, "Budget", "r0")
POLICIES = Feature("policies", int, 2, "Policies", "K")

# SOLVER is REQUIRED.
SOLVER = Feature("solver", Solver, MIP, "Solver", "Sol", SOLVERS)
# Solver Parameter Features - they define the PARAMETERS OF THE SOLVER (we can repeat identical runs, with different solver parameters).
# None of these are required.
SUBSOLVER = Feature("subsolver", Solver, NONE, "Sub Solver", "Sub", SOLVERS)
M_SYM = Feature("m_sym", int, -1, "Manual Symmetry", "Msym")
G_SYM = Feature("g_sym", int, -1, "Gurobi Symmetry", "Gsym")
# Instance Features + Solver + Solver Parameters = a Run.


# Output Features - they define the OUTPUT OF THE RUN.
OBJECTIVE = Feature("objective", float, 0.0, "Objective", "Obj")
UNBOUNDED = Feature("unbounded", str, "NOT_UNBOUNDED", "Unbounded", "Unb")
OPTIMAL = Feature("optimal", str, "OPTIMAL", "Optimal", "Opt")
PARTITION = Feature("partition", str, "0-1-1-0-0", "Partition", "Part")
CUTS_ROUNDS = Feature("cuts_rounds", int, 0, "Callbacks", "Cb")
CUTS_ADDED = Feature("cuts_added", int, 0, "Cuts Added", "Cuts")
GAP = Feature("gap", float, 0.0, "Gap (%)", "Gap (%)")
# TODO: could be a special time type.
AVG_CBTIME = Feature("avg_cbtime", float, 0.0, "Average Callback Time (ms)", "CbT (ms)")
AVG_SPTIME = Feature("avg_sptime", float, 0.0, "Average Cut Separation Time (ms)", "CutT (ms)")
TIME = Feature("time", float, 0.0, "Running Time (ms)", "T (ms)")
AVG_CBTIME_S = Feature("avg_cbtime_s", float, 0.0, "Average Callback Time (s)", "CbT (s)")
AVG_SPTIME_S = Feature("avg_sptime_s", float, 0.0, "Average Cut Separation Time (s)", "CutT (s)")
TIME_S = Feature("time_s", float, 0.0, "Running Time (s)", "T (s)")
# TODO: a special feature that has a default value of -1, but is a float.
EMPIRICAL_OPTIMAL_RATIO = Feature("empirical_optimal_ratio", float, -1.0, "Empirical Optimal Ratio", "rOpt")
EMPIRICAL_SUBOPTIMAL_RATIO = Feature("empirical_suboptimal_ratio", float, -1.0, "Empirical Suboptimal Ratio", "rObj")
BEST_OPTIMAL = Feature("best_optimal", float, -1.0, "Best Optimal Objective", "MaxOpt")
BEST_OBJECTIVE = Feature("best_objective", float, -1.0, "Best Objective", "MaxObj")
EXACT_ALPHA = Feature("exact_alpha", float, -1.0, "Exact Alpha", "a_exact")
EXACT_ALPHA_TIME_S = Feature("exact_alpha_time_s", float, -1.0, "Exact Alpha Time (s)", "a_e_T (s)")
ALPHA_HAT_ONE = Feature("alpha_hat_one", float, -1.0, "Alpha Hat One", "a_hat1")
ALPHA_HAT_ONE_TIME_S = Feature("alpha_hat_one_time_s", float, -1.0, "Alpha Hat One Time (s)", "a_hat1_T (s)")
ALPHA_HAT_TWO = Feature("alpha_hat_two", float, -1.0, "Alpha Hat Two", "a_hat2")
ALPHA_HAT_TWO_TIME_S = Feature("alpha_hat_two_time_s", float, -1.0, "Alpha Hat Two Time (s)", "a_hat2_T (s)")
UNINTERDICTED_SHORTEST_PATH = Feature("uninterdicted_shortest_path", float, -1.0, "Uninterdicted Shortest Path", "USP")
ADAPTIVE_INCREMENT = Feature("adaptive_increment", float, -1.0, "Adaptive Increment", "a_inc")
# Run + Outputs = Row.


MIP.add_parameters([M_SYM, G_SYM])
MIP.add_latex_output_features([])
BENDERS.add_parameters([M_SYM, G_SYM])
BENDERS.add_latex_output_features([])
ENUMERATION.add_parameters([SUBSOLVER])
ENUMERATION.add_latex_output_features([])
GREEDY.add_parameters([SUBSOLVER])
GREEDY.add_latex_output_features([])


################
##### TODO: STARTING HERE, STILL NEEDS TO BE REFACTORED

COLS = {
    "name_inputs_str": [
        SET_NAME,
        INSTANCE_NAME,
    ],
    "graph_inputs_int": [
        NODES,
        ARCS,
        K_ZERO,
    ],
    "graph_inputs_rat": [
        DENSITY,
    ],
    "cost_instance_inputs_int": [
        SCENARIOS,
    ],
    "run_input_params_int": [
        BUDGET,
        POLICIES,
    ],
    "run_input_hyperparams_int": [
        M_SYM,
        G_SYM,
    ],
    "run_input_hyperparams_cat": [
        SUBSOLVER,
    ],
    "run_input_params_cat": [
        SOLVER,
    ],
    "solution_outputs_rat": [
        OBJECTIVE,
    ],
    "solution_outputs_cat": [
        UNBOUNDED,
        OPTIMAL,
        PARTITION,
    ],
    "model_outputs_int": [
        CUTS_ROUNDS,
        CUTS_ADDED,
    ],
    "model_outputs_rat": [
        GAP,
    ],
    "time_outputs_rat": [
        AVG_CBTIME,
        AVG_SPTIME,
        TIME,
    ],
    "time_outputs_s_rat": [
        AVG_CBTIME_S,
        AVG_SPTIME_S,
        TIME_S,
    ],
    "slow_constants": [
        EMPIRICAL_OPTIMAL_RATIO,
        EMPIRICAL_SUBOPTIMAL_RATIO,
        BEST_OPTIMAL,
        BEST_OBJECTIVE,
        EXACT_ALPHA,
        EXACT_ALPHA_TIME_S,
        ALPHA_HAT_ONE,
        ALPHA_HAT_ONE_TIME_S,
        ALPHA_HAT_TWO,
        ALPHA_HAT_TWO_TIME_S,
        UNINTERDICTED_SHORTEST_PATH,
        ADAPTIVE_INCREMENT,
    ],
}

FEATURES = []
for group in COLS.values():
    for f in group:
        FEATURES.append(f)

COLS["name_inputs"] = COLS["name_inputs_str"]
COLS["graph_inputs"] = COLS["graph_inputs_int"] + COLS["graph_inputs_rat"]
COLS["cost_instance_inputs"] = COLS["cost_instance_inputs_int"]
COLS["run_input_params"] = COLS["run_input_params_int"] + COLS["run_input_params_cat"]
COLS["run_input_hyperparams"] = (
    COLS["run_input_hyperparams_int"] + COLS["run_input_hyperparams_cat"]
)
COLS["run_inputs"] = COLS["run_input_params"] + COLS["run_input_hyperparams"]

COLS["inputs"] = (
    COLS["name_inputs"]
    + COLS["graph_inputs"]
    + COLS["cost_instance_inputs"]
    + COLS["run_inputs"]
)

COLS["model_outputs"] = COLS["model_outputs_int"] + COLS["model_outputs_rat"]

COLS["outputs_int"] = COLS["model_outputs_int"]
COLS["outputs_rat"] = COLS["model_outputs_rat"] + COLS["time_outputs_rat"]
COLS["outputs_cat"] = COLS["solution_outputs_cat"]

COLS["solution_outputs"] = COLS["solution_outputs_rat"] + COLS["solution_outputs_cat"]
COLS["outputs"] = (
    COLS["solution_outputs"] + COLS["model_outputs"] + COLS["time_outputs_rat"]
)


COLS["raw"] = COLS["inputs"] + COLS["outputs"]
# "processed" by slow computations in add_ratios.py. (persistent / db)
COLS["processed"] = COLS["raw"] + COLS["slow_constants"]
# "finished" by fast clean up in common_cleanup() (local to callstack)
COLS["logging_run_header"] = COLS["name_inputs_str"] + COLS["run_inputs"]
COLS["logging_outputs"] = COLS["solution_outputs"] + COLS["time_outputs_s_rat"]
COLS["finished"] = COLS["processed"] + COLS["time_outputs_s_rat"]

# Additional categorizations (no new columns past this point).
# Grouping - column subsets useful for common groupby operations, such as finding dups.
COLS["same_instance"] = (
    COLS["name_inputs"] + COLS["graph_inputs"] + COLS["cost_instance_inputs"]
)
# TODO REORGANIZE THESE CONCEPTS
COLS["same_run"] = COLS["same_instance"] + COLS["run_input_params_int"]
COLS["same_run_params"] = COLS["same_instance"] + COLS["run_input_params"]
COLS["same_run_hyperparams"] = COLS["inputs"]
COLS["same_run_and_outputs"] = (
    COLS["inputs"] + COLS["model_outputs"] + COLS["time_outputs_rat"]
)

# Types - column subsets partitioned by type, useful for cleaning functions - including all columns in "finished" - these subsets should be used in common_cleanup() once we are at "finished" columns.
COLS["integer"] = (
    COLS["graph_inputs_int"]
    + COLS["cost_instance_inputs_int"]
    + COLS["run_input_params_int"]
    + COLS["run_input_hyperparams_int"]
)
COLS["rational"] = (
    COLS["graph_inputs_rat"]
    + COLS["outputs_rat"]
    + COLS["time_outputs_s_rat"]
    + COLS["slow_constants"]
)

DP = {key.name: 2 for key in COLS["rational"]}
#########
