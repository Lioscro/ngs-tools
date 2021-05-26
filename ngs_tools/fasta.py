import re
from typing import Dict

from . import sequence, utils
from .logging import logger

from typing_extensions import Literal


class FastaEntryError(Exception):
    pass


class FastaError(Exception):
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


def split_genomic_fasta_to_cdna(
    fasta_path: str, out_path: str, gene_infos: dict, transcript_infos: dict
) -> str:
    """Split a genomic FASTA into cDNA by using gene and transcript information
    generated from extracting information from a GTF.

    Args:
        fasta_path: Path to FASTA containing genomic sequences
        out_path: Path to output FASTA that will contain cDNA sequences
        gene_infos: Dictionary containing gene information, as returned by
            :func:`ngs_tools.gtf.genes_and_transcripts_from_gtf`
        transcript_infos: Dictionary containing transcript information, as returned by
            :func:`ngs_tools.gtf.genes_and_transcripts_from_gtf`

    Returns:
        Path to written FASTA
    """
    with Fasta(fasta_path, 'r') as f_in, Fasta(out_path, 'w') as f_out:
        for entry in f_in:
            # Find all gene and transcripts in this chromosome
            _gene_infos = {}
            _transcript_infos = {}
            for gene_id, gene_attributes in gene_infos.items():
                if gene_attributes['chromosome'] == entry.name:
                    _gene_infos[gene_id] = gene_attributes
                    _transcript_infos.update({
                        transcript_id: transcript_infos[transcript_id]
                        for transcript_id in gene_attributes['transcripts']
                    })

            # Write all transcripts as separate FASTA entries.
            for transcript_id, transcript_attributes in _transcript_infos.items(
            ):
                gene_id = transcript_attributes['gene_id']
                gene_name = _gene_infos[gene_id].get('gene_name')
                transcript_name = transcript_attributes.get('transcript_name')
                chromosome = _gene_infos[gene_id]['chromosome']
                segment = transcript_attributes['segment']
                strand = _gene_infos[gene_id]['strand']
                header = FastaEntry.make_header(
                    transcript_id, {
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'transcript_name': transcript_name,
                        'chr': chromosome,
                        'start': segment.start + 1,
                        'end': segment.end,
                        'strand': strand
                    }
                )

                # Construct sequence
                s = ''
                for exon in transcript_attributes['exons']:
                    s += entry.sequence[exon.start:exon.end]
                if s:
                    if strand == '-':
                        s = sequence.complement_sequence(s, reverse=True)
                    cdna_entry = FastaEntry(header, s)
                    f_out.write(cdna_entry)
    return out_path


def split_genomic_fasta_to_intron(
    fasta_path: str,
    out_path: str,
    gene_infos: dict,
    transcript_infos: dict,
    flank: int = 30
) -> str:
    """Split a genomic FASTA into introns by using gene and transcript information
    generated from extracting information from a GTF. Optionally append flanking
    sequences and collapse introns that have overlapping flanking regions.

    Args:
        fasta_path: Path to FASTA containing genomic sequences
        out_path: Path to output FASTA that will contain cDNA sequences
        gene_infos: Dictionary containing gene information, as returned by
            :func:`ngs_tools.gtf.genes_and_transcripts_from_gtf`
        transcript_infos: Dictionary containing transcript information, as returned by
            :func:`ngs_tools.gtf.genes_and_transcripts_from_gtf`
        flank: Number of flanking bases to include for each intron. Defaults to 30.

    Returns:
        Path to written FASTA
    """
    with Fasta(fasta_path, 'r') as f_in, Fasta(out_path, 'w') as f_out:
        for entry in f_in:
            # Find all gene and transcripts in this chromosome
            _gene_infos = {}
            _transcript_infos = {}
            for gene_id, gene_attributes in gene_infos.items():
                if gene_attributes['chromosome'] == entry.name:
                    _gene_infos[gene_id] = gene_attributes
                    _transcript_infos.update({
                        transcript_id: transcript_infos[transcript_id]
                        for transcript_id in gene_attributes['transcripts']
                    })

            # Write all transcripts as separate FASTA entries.
            for transcript_id, transcript_attributes in _transcript_infos.items(
            ):
                gene_id = transcript_attributes['gene_id']
                gene_name = _gene_infos[gene_id].get('gene_name')
                transcript_name = transcript_attributes.get('transcript_name')
                chromosome = _gene_infos[gene_id]['chromosome']
                segment = transcript_attributes['segment']
                strand = _gene_infos[gene_id]['strand']

                # Create new SegmentCollection with appropriate flanks
                flanked = transcript_attributes['introns'].flank(
                    flank, left=segment.start, right=segment.end
                )

                for i, intron in enumerate(flanked):
                    header = FastaEntry.make_header(
                        f'{transcript_id}-I.{i+1}', {
                            'gene_id': gene_id,
                            'gene_name': gene_name,
                            'transcript_name': transcript_name,
                            'chr': chromosome,
                            'start': segment.start + 1,
                            'end': segment.end,
                            'strand': strand
                        }
                    )
                    s = entry.sequence[intron.start:intron.end]
                    if strand == '-':
                        s = sequence.complement_sequence(s, reverse=True)
                    intron_entry = FastaEntry(header, s)
                    f_out.write(intron_entry)

    return out_path
