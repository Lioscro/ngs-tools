from abc import ABC
from typing import Any, List

from .ArgumentValidator import ArgumentValidator, NoValidator


class ArgumentError(Exception):
    pass


class Argument(ABC):
    """Abstract class representing a single argument.
    """

    def __init__(
        self,
        name: str,
        pre_validator: ArgumentValidator = NoValidator(),
        post_validator: ArgumentValidator = NoValidator(),
        required: bool = True,
    ):
        self._name = name
        self._pre_validator = pre_validator
        self._post_validator = post_validator
        self._required = required

    @property
    def name(self) -> str:
        return self._name

    @property
    def required(self) -> bool:
        return self._required

    def pre_execute(self, arg: str):
        if not self._pre_validator(arg):
            raise ArgumentError(
                f'Value `{arg}` for argument `{self.name}` failed pre-execute validation '
                f'using validator `{self._pre_validator}`.'
            )

    def post_execute(self, arg: str):
        if not self._post_validator(arg):
            raise ArgumentError(
                f'Value `{arg}` for argument `{self.name}` failed post-execute validation '
                f'using validator `{self._post_validator}`.'
            )

    def render(self, arg: str) -> List[str]:
        return [arg]


# This is just a regular argument.
class PositionalArgument(Argument):
    pass


class ConstantArgument(Argument):

    def render(self, arg: Any) -> List[str]:
        return [self.name]


class NamedArgument(Argument):

    def render(self, arg: str) -> List[str]:
        return [self.name, arg]
