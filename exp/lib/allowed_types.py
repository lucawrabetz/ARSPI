from typing import Type, Any, Dict
from optiface.lib.solver_type import PySolverType

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
            SolverType: 'PySolverType',
        }
        self.string_to_type: Dict[str, Type[Any]] = {
            'int': int,
            'float': float,
            'str': str,
            'bool': bool,
            'SolverType': PySolverType,
        }
        self.type_to_conversion: Dict[Type[Any], Callable[str, Any]] = {
            int: int,
            float: float,
            str: lambda x: x,
            bool: convert_bool(),
            PySolverType: convert_solver_type(),
        }

    def convert_bool(self, default: str) -> bool:
        if default.lower() == 'true':
            return True
        elif default.lower() == 'false':
            return False
        else:
            raise ValueError(f"Default value {default} is not a valid boolean")

    def convert_solver_type(self, default: str) -> PySolverType:
        # TODO this should be a simple access to a Dict of solver type string names to their declared PySolverType objects. This dict has to be declared by the IMPLEMENTER, close to where they declare the solver types features and feature types.
        pass

    def convert_default(self, default: str, data_type: Type[Any]) -> Any:
        """
        Convert the string default value to the correct type.
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

_ALLOWED_TYPES = TypeMapper()

