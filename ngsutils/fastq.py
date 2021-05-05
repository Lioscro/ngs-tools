import gzip
import os
from typing import Generator, Optional, Union

import pysam
import shortuuid
from tqdm import tqdm

from . import utils


class ReadError(Exception):
    pass


class FastqError(Exception):
    pass


class Quality:

    def __init__(self, qualities: str):
        self._string = qualities

    @property
    def string(self):
        return self._string

    @property
    def values(self):
        return list(pysam.qualitystring_to_array(self._string))

    @property
    def probs(self):
        return [10**(-q / 10) for q in self.values]


class Read:
    """Class that represents a FASTQ read. Once the class is initialized, with the
    header, sequence, and quality strings, they should not be changed (hence using @property).
    Only Phred33 quality scores are supported, but there is no check for this.
    """

    def __init__(self, header: str, sequence: str, qualities: str):
        if not header.startswith('@'):
            raise ReadError(f'FASTQ header `{header}` does not start with `@`')

        self._header = header.strip()
        self._sequence = sequence.strip()
        self._qualities = Quality(qualities.strip())

    @property
    def header(self):
        return self._header

    @property
    def sequence(self):
        return self._sequence

    @property
    def qualities(self):
        return self._qualities

    @property
    def name(self):
        return self.header[1:self.header.
                           index(' ')] if ' ' in self.header else self.header

    @property
    def extras(self):
        return self.header[self.header.index(' ') +
                           1:] if ' ' in self.header else ''


class Fastq:
    """Class that represents a FASTQ file. The file is opened on initialization
    and closed when the object is garbage collected (i.e. deleted).
    """

    def __init__(self, path: str, mode: str = 'r'):
        self.path = path
        self.mode = mode
        self.fp = None
        self.closed = False

        # Immediately open file descriptor
        self._open()

    def __del__(self):
        self.close()

    @property
    def is_gzip(self):
        if os.path.exists(self.path):
            return utils.is_gzip(self.path)
        else:
            return self.path.endswith('.gz')

    def _open(self):
        self.fp = (
            gzip.open(self.path, f'{self.mode}t')
            if self.is_gzip else open(self.path, self.mode)
        )

    def close(self):
        self.closed = True
        self.fp.close()

    def reset(self):
        """Reset the internal file descriptor to the beginning of the file.
        """
        self.fp.seek(0)

    def read(self) -> Union[Read, None]:
        """Read a single FASTQ entry.
        """
        if self.mode != 'r':
            raise FastqError(f'Can not read from file in mode `{self.mode}`')
        if self.closed:
            raise FastqError('Can not read from closed file')

        header, sequence, _, qualities = [next(self.fp) for _ in range(4)]
        return Read(header, sequence, qualities)

    def write(self, read: Read):
        """Write a single FASTQ entry.
        """
        if self.mode != 'w':
            raise FastqError(f'Can not write to file in mode `{self.mode}`')
        if self.closed:
            raise FastqError('Can not write to closed file')

        self.fp.write(f'{read.header}\n')
        self.fp.write(f'{read.sequence}\n')
        self.fp.write('+\n')
        self.fp.write(f'{read.qualities.string}\n')

    def reads(self, n: Optional[int] = None) -> Generator[Read, None, None]:
        """Generator to read n reads from the Fastq file as a list of Reads.
        If n is not provided, all reads are yielded (starting from the current
        location of the file descriptor).
        """
        while True:
            try:
                read = self.read()
            except StopIteration:
                return
            yield read


def fastq_to_bam(
    fastq_path: str,
    bam_path: str,
    name: Optional[str] = None,
    n_threads: int = 1
):
    """Convert a Fastq to unmapped BAM.
    """
    fastq = Fastq(fastq_path)
    rg = name or shortuuid.uuid()
    header = pysam.AlignmentHeader.from_dict({
        'HD': {
            'VN': pysam.version.__samtools_version__,
            'SO': 'unsorted'
        },
        'RG': [{
            'ID': rg
        }],
    })
    with pysam.AlignmentFile(bam_path, 'wb', header=header,
                             threads=n_threads) as f:
        for read in tqdm(fastq.reads(), smoothing=0, desc='Writing BAM'):
            al = pysam.AlignedSegment(header)
            al.query_name = read.name
            al.query_sequence = read.sequence
            al.query_qualities = read.qualities.values
            al.flag = 4  # unmapped
            al.tags = [('RG', rg)]
            f.write(al)
