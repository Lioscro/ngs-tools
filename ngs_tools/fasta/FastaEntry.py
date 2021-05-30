import re
from typing import Dict

from ..logging import logger


class FastaEntryError(Exception):
    pass


class FastaEntry:
    """Represents a single FASTA entry, which consists of a header and a sequence.

    Attributes:
        ATTRIBUTE_PARSER: Static attribute that is a compiled regex. Used to parse
            attributes.
        _header: Header string; for internal use only. Use :attr:`header` instead.
        _sequence: Sequence string; for internal use only. Use :attr:`sequence` instead.
    """
    ATTRIBUTE_PARSER = re.compile(r'(?P<key>\S+?):(?P<value>\S*)')

    def __init__(self, header: str, sequence: str):
        """
        Args:
            header: Header string, including the ``>`` character
            sequence: Sequence string

        Raises:
            FastaEntryError: if the ``header`` does not start with ``>``
        """
        if not header.startswith('>'):
            raise FastaEntryError(
                f'FASTA header `{header}` does not start with `>`'
            )

        self._header = header.strip()
        self._sequence = sequence.strip()

    @property
    def header(self) -> str:
        """Header string, including the ``>`` character"""
        return self._header

    @property
    def sequence(self) -> str:
        """Sequence string"""
        return self._sequence

    @property
    def name(self) -> str:
        """Name of the sequence, which comes immedately after ``>`` in the header
        """
        return self.header[1:self.header.index(' ')
                           ] if ' ' in self.header else self.header[1:]

    @property
    def attributes(self) -> Dict[str, str]:
        """Dictionary of entry attributes, parsed from the substring of the header
        after the first space character.
        """
        attributes = {}
        attribute_str = self.header[self.header.index(' ') +
                                    1:] if ' ' in self.header else ''
        for key, value in self.ATTRIBUTE_PARSER.findall(attribute_str):
            if key in attributes:
                logger.warning(
                    f'Duplicate key in FASTA entry {self.name}: {key}'
                )
            attributes[key] = value
        return attributes

    @staticmethod
    def make_header(name: str, attributes: Dict[str, str]) -> str:
        """Static method to construct a header string from a name and attributes.

        Args:
            name: entry name
            attributes: dictionary containing entry attributes
        """
        attributes_str = ' '.join(f'{k}:{v}' for k, v in attributes.items())
        return f'>{name} {attributes_str}'
