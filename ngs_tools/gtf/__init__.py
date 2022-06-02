from itertools import product
from typing import Callable, Dict, Generator, List, Optional, Set, Tuple, Union

from ..logging import logger
from ..progress import progress
from .Gtf import Gtf, GtfError
from .GtfEntry import GtfEntry, GtfEntryError
from .Segment import Segment, SegmentError
from .SegmentCollection import SegmentCollection, SegmentCollectionError


def parse_gtf(
    gtf_path: str,
    filter_func: Callable[[GtfEntry], bool] = lambda entry: True,
    show_progress: bool = False,
) -> Generator[GtfEntry, None, None]:
    """Parse GTF and yield only the specified features as :class:`GtfEntry` instances.

    Args:
        gtf_path: path to GTF file
        filter_func: Function that takes a :class:`GtfEntry` instance and
            returns ``True`` for entries to process and ``False`` for entries
            to ignore. Defaults to no filtering.
        show_progress: Whether to display a progress bar. Defaults to False.

    Yields:
        GTF entries
    """
    with Gtf(gtf_path, 'r') as f:
        for entry in progress(f, desc='Parsing GTF', disable=not show_progress):
            if not filter_func(entry):
                continue
            yield entry


def genes_and_transcripts_from_gtf(
    gtf_path: str,
    use_version: bool = False,
    filter_func: Callable[[GtfEntry], bool] = lambda entry: True,
    show_progress: bool = False,
) -> Tuple[dict, dict]:
    """Parse GTF for gene and transcript information. Also, compute the introns of
    each transcript.

    Args:
        gtf_path: path to GTF file
        use_version: whether or not to use gene and transcript versions
        filter_func: Function that takes a :class:`GtfEntry` instance and
            returns ``True`` for entries to process and ``False`` for entries
            to ignore. Defaults to no filtering.
        show_progress: Whether to display a progress bar. Defaults to False.

    Returns:
        Dictionary containing gene information
        Dictionary containing transcript information
    """
    gene_infos = {}
    transcript_exons = {}
    transcript_infos = {}
    renamed = set()

    entry_filter = lambda entry: (
        entry.feature in ['exon', 'transcript', 'gene']
    ) and filter_func(entry)

    for entry in parse_gtf(gtf_path, filter_func=entry_filter,
                           show_progress=show_progress):
        # IMPORTANT: every feature must have gene_id
        attributes = entry.attributes
        gene_id = attributes.get('gene_id')
        if not gene_id:
            logger.warning(
                f'Found feature `{entry.feature}` without `gene_id`. Ignoring.'
            )
            continue
        if use_version:
            gene_version = attributes.get('gene_version')
            if gene_version:
                gene_id = f'{gene_id}.{gene_version}'
        gene_name = attributes.get('gene_name', '')

        # Update gene info. If we have a gene feature, override instead
        formatted = {
            'segment': Segment(entry.start - 1, entry.end),
            'chromosome': entry.chromosome,
            'strand': entry.strand,
        }
        if entry.feature == 'gene':
            gene_infos.setdefault(gene_id, {}).update(formatted)
            if 'transcripts' not in gene_infos[gene_id]:
                gene_infos[gene_id]['transcripts'] = []
        gene_info = gene_infos.setdefault(gene_id, formatted)
        segment = gene_info['segment']
        gene_info['segment'] = Segment(
            min(segment.start, entry.start - 1), max(segment.end, entry.end)
        )
        gene_info['gene_name'] = gene_name or gene_info.get('gene_name', '')
        gene_info.setdefault('transcripts', [])

        # Transcript and exon features
        # Some GTFs (specifically gencode) add transcript_id to gene features,
        # but in these cases gene_id = transcript_id, so we ignore these explicitly.
        if entry.feature != 'gene':
            # IMPORTANT: every transcript, exon feature must have transcript_id
            transcript_id = attributes.get('transcript_id')
            if not transcript_id:
                logger.warning(
                    f'Gene `{gene_id}` feature `{entry.feature}` does not have a transcript_id. Ignoring.'
                )
                continue
            transcript_name = attributes.get('transcript_name', '')
            if use_version:
                transcript_version = attributes.get('transcript_version')
                if transcript_version:
                    transcript_id = f'{transcript_id}.{transcript_version}'

            # Check for duplicates and rename appropriately
            # The assumption is gene_id, transcript_id pairs MUST be unique
            # (i.e. there are no distinct transcripts with the same ID within a gene)
            if gene_id != transcript_infos.get(transcript_id, {}).get('gene_id',
                                                                      gene_id):
                logger.warning(
                    f'Transcript `{transcript_id}` is assigned to multiple genes. '
                    'All instances of this transcript will be renamed in the following format: '
                    f'`{gene_id}-{transcript_id}`.'
                )
                renamed.add(transcript_id)
                _transcript_id = f'{gene_id}-{transcript_id}'
                # Retroactively rename previously parsed transcript_id
                # Note that any transcript IDs that don't have an entry in transcript_infos
                # will be removed in the next stage.
                gene_infos[transcript_infos[transcript_id]['gene_id']
                           ]['transcripts'].append(_transcript_id)
                transcript_infos[_transcript_id] = transcript_infos[
                    transcript_id].copy()
                del transcript_infos[transcript_id]
                if transcript_id in transcript_exons:
                    transcript_exons[_transcript_id] = transcript_exons[
                        transcript_id].copy()
                    del transcript_exons[transcript_id]

            if transcript_id in renamed:
                transcript_id = f'{gene_id}-{transcript_id}'

            gene_info['transcripts'].append(transcript_id)
            formatted = {
                'gene_id': gene_id,
                'segment': Segment(entry.start - 1, entry.end),
            }
            if entry.feature == 'transcript':
                transcript_infos.setdefault(transcript_id, {}).update(formatted)
            transcript_info = transcript_infos.setdefault(
                transcript_id, formatted
            )
            segment = transcript_info['segment']
            transcript_info['segment'] = Segment(
                min(segment.start, entry.start - 1),
                max(segment.end, entry.end)
            )
            transcript_info['transcript_name'
                            ] = transcript_name or transcript_info.get(
                                'transcript_name', ''
                            )

            if entry.feature == 'exon':
                transcript_exons.setdefault(
                    transcript_id, SegmentCollection()
                ).add_segment(Segment(entry.start - 1, entry.end))

    # Clean gene infos so that they link to only transcripts existing in transcript_infos
    for gene_id, attributes in gene_infos.items():
        cleaned = list(
            set(t for t in attributes['transcripts'] if t in transcript_infos)
        )
        attributes['transcripts'] = cleaned
        if not cleaned:
            logger.warning(
                f'Gene `{gene_id}` has no transcripts. '
                f'The entire gene will be marked as a transcript and an exon with ID `{gene_id}`.'
            )
            attributes['transcripts'].append(gene_id)
            transcript_infos[gene_id] = {
                'gene_id': gene_id,
                'transcript_name': '',
                'segment': attributes['segment'],
            }
            transcript_exons[gene_id] = SegmentCollection(
                segments=[attributes['segment']]
            )

    # Calculate introns
    for transcript_id, attributes in transcript_infos.items():
        transcript_interval = attributes['segment']

        exons = transcript_exons.get(transcript_id, SegmentCollection())
        introns = exons.invert(transcript_interval)
        if not exons:
            logger.warning(
                f'Gene `{attributes["gene_id"]}` transcript `{transcript_id}` has no exons. '
                'The entire transcript will be marked as an intron.'
            )

        attributes['exons'] = exons
        attributes['introns'] = introns

    logger.debug(
        f'Parsed {len(gene_infos)} genes and {len(transcript_infos)} transcripts'
    )
    return gene_infos, transcript_infos
