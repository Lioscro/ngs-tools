from .. import utils
from .GtfEntry import GtfEntry


class GtfError(Exception):
    pass


class Gtf(utils.FileWrapper):
    """Class that represents a GTF file.
    """

    def read(self) -> GtfEntry:
        """Read a single GTF entry as a :class:`GtfEntry` instance.

        Returns:
            The next GTF entry

        Raises:
            GtfError: If the file was not opened for reading, or the file was closed.
            StopIteration: When there are no more entries to read.
        """
        if self.mode != 'r':
            raise GtfError(f'Can not read from file in mode `{self.mode}`')
        if self.closed:
            raise GtfError('Can not read from closed file')

        while True:
            line = next(self.fp).strip()
            if not line.startswith('#'):
                break

        return GtfEntry(line)

    def write(self, entry: GtfEntry):
        """Write a single GTF entry.

        Args:
            entry: The GTF entry to write

        Raises:
            FastaError: If the file was not opened for writing, or the file was closed.
        """
        if self.mode != 'w':
            raise GtfError(f'Can not write to file in mode `{self.mode}`')
        if self.closed:
            raise GtfError('Can not write to closed file')

        self.fp.write(f'{entry.line}\n')
