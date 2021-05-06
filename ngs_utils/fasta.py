import re
from typing import Literal

from . import utils
from .logging import logger


class FastaEntryError(Exception):
    pass


class FastaError(Exception):
    pass


class FastaEntry:
    ATTRIBUTE_PARSER = re.compile(r'(?P<key>\S+?):(?P<value>\S*)')

    def __init__(self, header: str, sequence: str):
        if not header.startswith('>'):
            raise FastaEntryError(
                f'FASTA header `{header}` does not start with `>`'
            )

        self._header = header.strip()
        self._sequence = sequence.strip()

    @property
    def header(self):
        return self._header

    @property
    def sequence(self):
        return self._sequence

    @property
    def name(self):
        return self.header[1:self.header.
                           index(' ')] if ' ' in self.header else self.header

    @property
    def attributes(self):
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
    def make_header(name, attributes):
        attributes_str = ' '.join(f'{k}:{v}' for k, v in attributes.items())
        return f'>{name} {attributes_str}'


class Fasta(utils.FileWrapper):
    PARSER = re.compile(r'^>(?P<sequence_id>\S+)(?P<group>.*)')
    GROUP_PARSER = re.compile(r'(?P<key>\S+?):(?P<value>\S+)')

    def __init__(self, path: str, mode: Literal['r', 'w'] = 'r'):
        super(Fasta, self).__init__(path, mode)

        # Cache for next header is needed to implement read() properly.
        self._header = None

    def read(self) -> FastaEntry:
        if self.mode != 'r':
            raise FastaError(f'Can not read from file in mode `{self.mode}`')
        if self.closed:
            raise FastaError('Can not read from closed file')

        # Read header
        header = self._header or next(self.fp)

        # Read rest as sequences until we find another > character
        sequence = ''
        while True:
            try:
                line = next(self.fp)
            except StopIteration:
                self._header = None
                break
            if line.startswith('>'):
                self._header = line
                break
            sequence += line

        return FastaEntry(header, sequence)

    def write(self, entry: FastaEntry):
        if self.mode != 'w':
            raise FastaError(f'Can not write to file in mode `{self.mode}`')
        if self.closed:
            raise FastaError('Can not write to closed file')

        self.fp.write(f'{entry.header}\n')
        self.fp.write(f'{entry.sequence}\n')
