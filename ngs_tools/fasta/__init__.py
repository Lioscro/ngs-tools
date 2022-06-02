from .. import sequence
from ..progress import progress
from .Fasta import Fasta, FastaError
from .FastaEntry import FastaEntry, FastaEntryError


def split_genomic_fasta_to_cdna(
    fasta_path: str,
    out_path: str,
    gene_infos: dict,
    transcript_infos: dict,
    show_progress: bool = False,
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
        show_progress: Whether to display a progress bar. Defaults to False.

    Returns:
        Path to written FASTA
    """
    with Fasta(fasta_path, 'r') as f_in, Fasta(out_path, 'w') as f_out:
        for entry in progress(f_in, desc='Splitting cDNA',
                              disable=not show_progress):
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
    flank: int = 30,
    show_progress: bool = False,
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
        show_progress: Whether to display a progress bar. Defaults to False.

    Returns:
        Path to written FASTA
    """
    with Fasta(fasta_path, 'r') as f_in, Fasta(out_path, 'w') as f_out:
        for entry in progress(f_in, desc='Splitting introns',
                              disable=not show_progress):
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
