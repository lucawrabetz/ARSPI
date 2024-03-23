from typing import final

from feature_interface import *
from util import *


RAW_FEATURES = {
    "name_inputs_str": [
        RawFeature("set_name", str, "layer"),
        RawFeature("instance_name", str, "layer_123"),
    ],
    "graph_inputs_int": [
        RawFeature("nodes", int, 123),
        RawFeature("arcs", int, 1231),
        RawFeature("k_zero", int, 5),
    ],
    "graph_inputs_rat": [
        RawFeature("density", float, 1231 / (123 * 122)),
    ],
    "cost_instance_inputs_int": [
        RawFeature("scenarios", int, 5),
    ],
    "run_input_params_int": [
        RawFeature("budget", int, 5),
        RawFeature("policies", int, 2),
    ],
    "run_input_hyperparams_int": [
        RawFeature("m_sym", int, -1),
        RawFeature("g_sym", int, -1),
    ],
    "run_input_hyperparams_cat": [
        RawFeature("subsolver", str, "NONE"),
    ],
    "run_input_params_cat": [
        RawFeature("solver", str, "MIP"),
    ],
    "solution_outputs_rat": [
        RawFeature("objective", float, 0.0),
    ],
    "solution_outputs_cat": [
        RawFeature("unbounded", str, "NOT_UNBOUNDED"),
        RawFeature("optimal", str, "OPTIMAL"),
        RawFeature("partition", str, "0-1-1-0-0"),
    ],
    "model_outputs_int": [
        RawFeature("cuts_rounds", int, 0),
        RawFeature("cuts_added", int, 0),
    ],
    "model_outputs_rat": [
        RawFeature("gap", int, 0.0),
    ],
    "time_outputs_rat": [
        RawFeature("avg_cbtime", int, 0.0),
        RawFeature("avg_sptime", int, 0.0),
        RawFeature("time", int, 0.0),
    ],
}
