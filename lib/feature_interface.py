from typing import Any
from abc import ABC, abstractmethod
from util import *


# Level 1
class Feature(ABC):
    DTYPE: type
    DEFAULT: Any

    def __init__(self):
        self.name = to_snakecase(self.__name__)


class SetName(Feature):
    pass


class InstanceName(Feature):
    pass


class Nodes(Feature):
    pass


class Arcs(Feature):
    pass


class KZero(Feature):
    pass


class Density(Feature):
    pass


class Scenarios(Feature):
    pass


class Budget(Feature):
    pass


class Policies(Feature):
    pass


class Solver(Feature):
    pass


class MSym(Feature):
    pass


class GSym(Feature):
    pass


class Subsolver(Feature):
    pass


class Objective(Feature):
    pass


class Unbounded(Feature):
    pass


class Optimal(Feature):
    pass


class Partition(Feature):
    pass


class CutsRounds(Feature):
    pass


class CutsAdded(Feature):
    pass


class Gap(Feature):
    pass


class AvgCbtime(Feature):
    pass


class AvgSptime(Feature):
    pass


class Time(Feature):
    pass


class EmpiricalOptimalRatio(Feature):
    pass


class EmpiricalSuboptimalRatio(Feature):
    pass


class BestOptimal(Feature):
    pass


class BestObjective(Feature):
    pass


class AvgCbtime_s(Feature):
    pass


class AvgSptime_s(Feature):
    pass


class Time_s(Feature):
    pass
