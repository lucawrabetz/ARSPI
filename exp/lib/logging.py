import warnings
# TODO: should this debugger class be a subclass of something bigger, or related to another logging class? See TODO in query.py asking for a logging class.
# Small debugger class that holds a debugger flag, and prints messages that are commonly used for debugging.
class Debugger:
    def __init__(self, debug: bool = False) -> None:
        self.debug = debug

    def log(self, message: str) -> None:
        if self.debug:
            print(message)

    def warn(self, message: str) -> None:
        if self.debug:
            warnings.warn(message)

DEBUG = Debugger(False)
