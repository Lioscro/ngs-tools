from typing_extensions import Literal

from .. import utils
from .FastaEntry import FastaEntry


class FastaError(Exception):
    pass


class Fasta(utils.FileWrapper):
    """Represents a single FASTA file.

    Attributes:
        _header: Variable that temporarily holds the header string for the next
            FASTA entry; for internal use only.
    """

    def __init__(self, path: str, mode: Literal['r', 'w'] = 'r'):
        super(Fasta, self).__init__(path, mode)

        # Cache for next header is needed to implement read() properly.
        self._header = None

    def read(self) -> FastaEntry:
        """Read a single FASTA entry as a :class:`FastaEntry` instance.

        Returns:
            The next FASTA entry

        Raises:
            FastaError: If the file was not opened for reading, or the file was closed.
            StopIteration: When there are no more entries to read.
        """
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
                line = next(self.fp).strip()
            except StopIteration:
                self._header = None
                break
            if line.startswith('>'):
                self._header = line
                break
            sequence += line

        return FastaEntry(header, sequence)

    def write(self, entry: FastaEntry):
        """Write a single FASTA entry.

        Args:
            entry: The FASTA entry to write

        Raises:
            FastaError: If the file was not opened for writing, or the file was closed.
        """
        if self.mode != 'w':
            raise FastaError(f'Can not write to file in mode `{self.mode}`')
        if self.closed:
            raise FastaError('Can not write to closed file')

        self.fp.write(f'{entry.header}\n')
        self.fp.write(f'{entry.sequence}\n')
