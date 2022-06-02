import os
import shutil
import subprocess
from typing import Dict, List, NamedTuple, Optional

from ..logging import logger
from ..utils import run_executable, set_executable


class BinaryError(Exception):
    pass


class BinaryResult(NamedTuple):
    command: List[str]
    process: subprocess.Popen
    stdout: Optional[str] = None
    stderr: Optional[str] = None


class Binary:
    """Wrapper around a binary file that provides a object-oriented interface to
    :func:`run_executable`.
    """

    def __init__(self, path: str):
        path = shutil.which(path)
        if not os.path.exists(path):
            raise BinaryError(f"{path} does not exist")
        path = os.path.abspath(path)
        if not os.access(path, os.X_OK):
            logger.warning(
                f"{path} is not executable. Attempting to change permissions."
            )
            set_executable(path)
        self._path = path

    @property
    def path(self) -> str:
        return self._path

    def __str__(self) -> str:
        return self.path

    def pre_execute(self, command: str, *args, **kwargs):
        pass

    def post_execute(self, command: str, result: BinaryResult, *args, **kwargs):
        pass

    def __call__(self, command: List[str], *args, **kwargs) -> BinaryResult:
        """Wrapper around the :func:`run_executable` function. The first element
        of the `command` argument to the function will always be the path to
        the current binary.
        """
        command = [self.path] + command
        self.pre_execute(command, *args, **kwargs)
        result = run_executable(command, *args, **kwargs)
        result = BinaryResult(command, *result) if isinstance(
            result, tuple
        ) else BinaryResult(command, result)
        self.post_execute(command, result, *args, **kwargs)
        return result


class BinaryExecutor:
    """Class that handles execution of binaries with specified arguments.
    """

    def __init__(self, binary: Binary, *args):
        self._binary = binary
        self._arguments = {}
        for arg in args:
            if arg.name in self._arguments:
                raise BinaryError(
                    f"Found duplicate argument named `{arg.name}`."
                )
            self._arguments[arg.name] = arg

    def pre_execute(self, command: str):
        pass

    def post_execute(self, command: str, result: BinaryResult):
        pass

    def __call__(self, values: Dict[str, Optional[str]], *args, **kwargs):
        """Run the binary with the specified values for each of the arguments.
        """
        arguments = self._arguments.copy()  # We will be popping each argument.
        command = []
        for name, value in values.items():
            if name not in arguments:
                raise BinaryError(
                    f'Unknown argument `{name}` for binary `{self._binary}`.'
                )
            arg = arguments.pop(name)
            arg.pre_execute(value)
            command.extend(arg.render(value))

        # Go through remaining and raise exception if any of them were required.
        for name, arg in arguments.items():
            if arg.required:
                raise BinaryError(
                    f'Argument `{arg.name}` is required for binary `{self._binary}`.'
                )

        # Execute!
        self.pre_execute(command)
        result = self._binary(command)
        self.post_execute(command, result)

        for name, value in values.items():
            arg = self._arguments[name]
            arg.post_execute(value)

        return result
