from typing import Optional

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
                           index(' ')] if ' ' in self.header else self.header[1:]

    @property
    def attributes(self):
        return self.header[self.header.index(' ') +
                           1:] if ' ' in self.header else ''


class Fastq(utils.FileWrapper):
    """Class that represents a FASTQ file. The file is opened on initialization
    and closed when the object is garbage collected (i.e. deleted).
    """

    def read(self) -> Read:
        """Read a single FASTQ entry.
        """
        if self.mode != 'r':
            raise FastqError(f'Can not read from file in mode `{self.mode}`')
        if self.closed:
            raise FastqError('Can not read from closed file')

        header, sequence, _, qualities = [next(self.fp).strip() for _ in range(4)]
        return Read(header, sequence, qualities)

    def write(self, entry: Read):
        """Write a single FASTQ entry.
        """
        if self.mode != 'w':
            raise FastqError(f'Can not write to file in mode `{self.mode}`')
        if self.closed:
            raise FastqError('Can not write to closed file')

        self.fp.write(f'{entry.header}\n')
        self.fp.write(f'{entry.sequence}\n')
        self.fp.write('+\n')
        self.fp.write(f'{entry.qualities.string}\n')


def fastq_to_bam(
    fastq_path: str,
    bam_path: str,
    name: Optional[str] = None,
    n_threads: int = 1
):
    """Convert a Fastq to unmapped BAM.
    """
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
    with Fastq(fastq_path,
               'r') as f_in, pysam.AlignmentFile(bam_path, 'wb', header=header,
                                                 threads=n_threads) as f_out:
        for read in tqdm(f_in, smoothing=0, desc='Writing BAM'):
            al = pysam.AlignedSegment(header)
            al.query_name = read.name
            al.query_sequence = read.sequence
            al.query_qualities = read.qualities.values
            al.flag = 4  # unmapped
            al.tags = [('RG', rg)]
            f_out.write(al)
