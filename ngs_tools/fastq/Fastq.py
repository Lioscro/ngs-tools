from .. import utils
from .Read import Read


class FastqError(Exception):
    pass


class Fastq(utils.FileWrapper):
    """Class that represents a FASTQ file.
    """

    def read(self) -> Read:
        """Read a single FASTQ entry.

        Returns:
            The next read as a :class:`Read` instance

        Raises:
            FastqError: If the file was not opened for reading, or the file was closed.
            StopIteration: When there are no more entries to read.
        """
        if self.mode != 'r':
            raise FastqError(f'Can not read from file in mode `{self.mode}`')
        if self.closed:
            raise FastqError('Can not read from closed file')

        header, sequence, _, qualities = [
            next(self.fp).strip() for _ in range(4)
        ]
        return Read(header, sequence, qualities)

    def write(self, entry: Read):
        """Write a single FASTQ entry.

        Args:
            entry: :class:`Read` instance to write to the FASTQ

        Raises:
            FastaError: If the file was not opened for writing, or the file was closed.
        """
        if self.mode != 'w':
            raise FastqError(f'Can not write to file in mode `{self.mode}`')
        if self.closed:
            raise FastqError('Can not write to closed file')

        self.fp.write(f'{entry.header}\n')
        self.fp.write(f'{entry.sequence}\n')
        self.fp.write('+\n')
        self.fp.write(f'{entry.qualities.string}\n')
