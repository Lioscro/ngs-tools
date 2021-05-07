import re
from typing import Literal

from . import sequence, utils
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
                           index(' ')] if ' ' in self.header else self.header[1:]

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
            for transcript_id, transcript_attributes in _transcript_infos.items():
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
                        'start': segment.start+1,
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
    fasta_path: str, out_path: str, gene_infos: dict, transcript_infos: dict, flank: int = 30
) -> str:
    """Split a genomic FASTA into introns by using gene and transcript information
    generated from extracting information from a GTF.
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
            for transcript_id, transcript_attributes in _transcript_infos.items():
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
                            'start': segment.start+1,
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
