import os
from abc import ABC, abstractmethod


class ArgumentValidator(ABC):
    """Abstract argument validator.
    """

    @abstractmethod
    def __call__(self, arg) -> bool:
        return False

    def __and__(self, other: 'ArgumentValidator') -> 'ArgumentValidator':
        return AndValidator(self, other)

    def __or__(self, other: 'ArgumentValidator') -> 'ArgumentValidator':
        return OrValidator(self, other)

    def __invert__(self) -> 'ArgumentValidator':
        return NotValidator(self)

    def __str__(self) -> str:
        return 'False'


class AndValidator(ArgumentValidator):

    def __init__(self, left: ArgumentValidator, right: ArgumentValidator):
        self._left = left
        self._right = right

    def __call__(self, arg) -> bool:
        return self._left(arg) and self._right(arg)

    def __str__(self) -> str:
        return f'({self._left} AND {self._right})'


class OrValidator(ArgumentValidator):

    def __init__(self, left: ArgumentValidator, right: ArgumentValidator):
        self._left = left
        self._right = right

    def __call__(self, arg) -> bool:
        return self._left(arg) or self._right(arg)

    def __str__(self) -> str:
        return f'({self._left} OR {self._right})'


class NotValidator(ArgumentValidator):

    def __init__(self, validator: ArgumentValidator):
        self._validator = validator

    def __call__(self, arg) -> bool:
        return not self._validator(arg)

    def __str__(self) -> str:
        return f'(NOT {self._validator})'


class NoValidator(ArgumentValidator):

    def __call__(self, arg) -> bool:
        return True

    def __str__(self) -> str:
        return 'True'


class IntegerValidator(ArgumentValidator):

    def __call__(self, arg) -> bool:
        try:
            v = int(arg)
            return v == float(v)
        except ValueError:
            return False

    def __str__(self) -> str:
        return 'IsInteger'


class FloatValidator(ArgumentValidator):

    def __call__(self, arg) -> bool:
        try:
            v = float(arg)
            return v != int(v)
        except ValueError:
            return False

    def __str__(self) -> str:
        return 'IsFloat'


class PositiveValidator(ArgumentValidator):

    def __call__(self, arg) -> bool:
        try:
            return float(arg) > 0
        except ValueError:
            return False

    def __str__(self) -> str:
        return 'IsPositive'


class FileValidator(ArgumentValidator):

    def __call__(self, arg) -> bool:
        return os.path.isfile(arg)

    def __str__(self) -> str:
        return 'IsFile'


class DirValidator(ArgumentValidator):

    def __call__(self, arg) -> bool:
        return os.path.isdir(arg)

    def __str__(self) -> str:
        return 'IsDir'


class ConstantValidator(ArgumentValidator):

    def __init__(self, constant: str):
        self._constant = constant

    def __call__(self, arg) -> bool:
        try:
            v = type(self._constant)(arg)
            return v == self._constant
        except ValueError:
            return False

    def __str__(self) -> str:
        return f'IsConstant({self._constant})'


IsInteger = IntegerValidator()
IsFloat = FloatValidator()
IsPositive = PositiveValidator()
IsPositiveInteger = IsInteger & IsPositive

IsDir = DirValidator()
IsFile = FileValidator()
