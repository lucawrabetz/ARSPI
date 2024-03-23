from util import *
from feature_interface import *


class DataModel(ABC):
    FEATURES: list[Feature] = []


class AspiDataModel(DataModel):
    FEATURES = [
        SetName(),
        InstanceName(),
        Nodes(),
        Arcs(),
    ]


# Leaving off:
# We are missing some connecting class or object to represent classifications or groups for the purpose of creating lists of column names to group on.
# So when we use this, one example would be to have a bool flag on every feature class:
# cols = [c.name for c in DataModel.FEATURES if c.flag == true]
# for _, group in df.groupby(cols):
#     pass

# Or we just use a list of just names for every group we want thats not a type, for now, then we can try and better abstract out certain groups, or all of them.
