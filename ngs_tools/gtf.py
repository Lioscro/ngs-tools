import bisect
import re
from itertools import product
from typing import Dict, Generator, List, Optional, Set, Tuple, Union

from . import utils
from .logging import logger


class SegmentError(Exception):
    pass


class Segment:
    """Class to represent an integer interval segment, zero-indexed, start-inclusive
    and end-exclusive.

    Attributes:
        _start: Segment start; for internal use only. Use :attr:`start` instead.
        _end: Segment end; for internal use only. Use :attr:`end` instead.
    """

    def __init__(self, start: int, end: int):
        """
        Args:
            start: Segment start
            end: Segment end

        Raises:
            SegmentError: If ``end <= start`` or ``start < 0``
        """
        if end <= start or start < 0:
            raise SegmentError(f'Invalid segment [{start}:{end})')

        self._start = start
        self._end = end

    @property
    def start(self) -> int:
        """Segment start"""
        return self._start

    @property
    def end(self) -> int:
        """Segment end"""
        return self._end

    @property
    def width(self) -> int:
        """Segment width"""
        return self.end - self.start

    def is_in(self, i: int) -> bool:
        """Evaluate whether an integer is contained within the segment.

        Args:
            i: Integer number to check

        Returns:
            True or False
        """
        return i >= self.start and i < self.end

    def is_exclusive(self, segment: 'Segment') -> bool:
        """Evaluate whether this segment is exclusive of another segment.

        Args:
            segment: :class:`Segment` object to check

        Returns:
            True or False
        """
        return self.end <= segment.start or self.start >= segment.end

    def is_overlapping(self, segment: 'Segment') -> bool:
        """Evaluate whether this segment overlaps with another segment.

        Args:
            segment: :class:`Segment` object to check

        Returns:
            True or False
        """
        return not self.is_exclusive(segment)

    def is_subset(self, segment: 'Segment') -> bool:
        """Evaluate whether this segment is a subset of another segment.

        Args:
            segment: :class:`Segment` object to check

        Returns:
            True or False
        """
        return self.start >= segment.start and self.end <= segment.end

    def is_superset(self, segment: 'Segment') -> bool:
        """Evaluate whether this segment is a superset of another segment.

        Args:
            segment: :class:`Segment` object to check

        Returns:
            True or False
        """
        return self.start <= segment.start and self.end >= segment.end

    def flank(
        self,
        l: int,  # noqa: E741
        left: Optional[int] = None,
        right: Optional[int] = None
    ) -> 'Segment':
        """Construct a new segment with start and end flanks of length ``l``.
        Optionally, limit the size of the new segment.

        Args:
            l: Flank length
            left: Clip the resulting segment's start to this value. Defaults to None.
                If not provided, no clipping is performed.
            right: Clip the resulting segment's end to this value. Defaults to None.
                If not provided, no clipping is performed.

        Returns:
            The new segment
        """
        left = left or 0
        right = right or self.end + l
        return Segment(max(left, self.start - l), min(self.end + l, right))

    def __iter__(self):
        return iter((self.start, self.end))

    def __eq__(self, other):
        return (self.start, self.end) == (other.start, other.end)

    def __lt__(self, other):
        return (self.start, self.end) < (other.start, other.end)

    def __gt__(self, other):
        return (self.start, self.end) > (other.start, other.end)

    def __repr__(self):
        return f'Segment {(self.start, self.end)}'


class SegmentCollection:
    """Class to represent a collection of integer interval segments, zero-indexed.
    The segments are always sorted and overlaps are collapsed.

    Attributes:
        segments: List of :class:`Segment` instances
    """

    def __init__(self, segments: Optional[List[Segment]] = None):
        """
        Args:
            segments: List of segments. Defaults to None.
        """
        self.segments = sorted(segments) if segments else []
        if self.segments:
            self.collapse()

    @property
    def start(self) -> int:
        """Leftmost value of all segments. 0 if there are no segments."""
        return self.segments[0].start if self.segments else 0

    @property
    def end(self) -> int:
        """Rightmost value of all segments. 0 if there are no segments."""
        return self.segments[-1].end if self.segments else 0

    def add_segment(self, segment: Segment):
        """Add a segment to the collection.

        Args:
            segment: Segment to add
        """
        bisect.insort_left(self.segments, segment)
        self.collapse()

    def add_collection(self, collection: 'SegmentCollection'):
        """Add a collection of segments to the collection.

        Args:
            collection: Collection to add
        """
        self.segments = sorted(self.segments + collection.segments)
        self.collapse()

    def __iter__(self):
        return iter(self.segments)

    def __getitem__(self, i):
        return self.segments[i]

    def __len__(self):
        return len(self.segments)

    def __bool__(self):
        return bool(self.segments)

    def collapse(self):
        """Collapse the segments in this collection such that there are no overlapping
        segments. Any overlapping segments are merged into a single large segment.
        """
        segments = []
        combined = None
        for segment in self.segments:
            if combined is None:
                combined = segment
                continue
            if combined.is_overlapping(segment):
                combined = Segment(
                    min(combined.start, segment.start),
                    max(combined.end, segment.end)
                )
            else:
                segments.append(combined)
                combined = segment
        if combined is not None:
            segments.append(combined)

        self.segments = segments

    def span_is_exclusive(self, collection: 'SegmentCollection') -> bool:
        """Evaluate whether the span of this collection is exclusive of that of
        another collection.

        Args:
            collection: :class:`SegmentCollection` object to check

        Returns:
            True or False
        """
        return self.end <= collection.start or self.start >= collection.end

    def is_overlapping(self, collection: 'SegmentCollection') -> bool:
        """Evaluate whether this collection overlaps with another collection.

        Args:
            collection: :class:`SegmentCollection` object to check

        Returns:
            True or False
        """
        if self.span_is_exclusive(collection):
            return False
        return any(
            segment_1.is_overlapping(segment_2)
            for segment_1, segment_2 in product(self, collection)
        )

    def is_subset(self, collection: 'SegmentCollection') -> bool:
        """Evaluate whether this collection is a subset of another collection.

        Args:
            collection: :class:`SegmentCollection` object to check

        Returns:
            True or False
        """
        if self.span_is_exclusive(collection):
            return False
        # Assume segments are sorted
        collection_i = 0
        for cls_segment in self.segments:
            found = False
            for collection_segment in collection.segments[collection_i:]:
                if cls_segment.is_subset(collection_segment):
                    found = True
                    break
                collection_i += 1

            if not found:
                return False
        return True

    def is_superset(self, collection: 'SegmentCollection') -> bool:
        """Evaluate whether this collection is a superset of another collection.

        Args:
            collection: :class:`SegmentCollection` object to check

        Returns:
            True or False
        """
        return collection.is_subset(self)

    def flank(
        self,
        l: int,  # noqa: E741
        left: Optional[int] = None,
        right: Optional[int] = None
    ) -> 'SegmentCollection':
        """Construct a new segment collection where all the segments have start and
        end flanks of length ``l``. Optionally, limit the span of the new segment collection.
        This is done by calling :func:`Segment.flank` on all the segments and initializing
        a new :class:`SegmentCollection`. Any overlaps are collapsed.

        Args:
            l: Flank length
            left: Clip the resulting collection's start to this value. Defaults to None.
                If not provided, no clipping is performed.
            right: Clip the resulting collection's end to this value. Defaults to None.
                If not provided, no clipping is performed.

        Returns:
            The new collection
        """
        segments = [segment.flank(l, left, right) for segment in self.segments]
        return SegmentCollection(segments)

    @classmethod
    def from_positions(
        cls, positions: Union[List[int], Set[int]]
    ) -> 'SegmentCollection':
        """Initialize a new collection given a list or set of integer positions.

        Args:
            positions: Integer positions to construct the collection from

        Returns:
            The new collection
        """
        collection = cls()
        positions = sorted(positions)
        start = positions[0]
        end = positions[0] + 1
        for i in range(1, len(positions)):
            left = positions[i - 1]
            right = positions[i]

            if right - left > 1:
                collection.add_segment(Segment(start, end))
                start = right
            end = right + 1

        collection.add_segment(Segment(start, end))
        return collection

    @classmethod
    def from_collections(
        cls, *collections: 'SegmentCollection'
    ) -> 'SegmentCollection':
        """Initialize a new collection given an arbitrary number of collections.

        Args:
            *collections: The collections to combine

        Returns:
            The new collection
        """
        segments = []
        for collection in collections:
            segments.extend(collection.segments)
        return cls(segments=segments)

    def __repr__(self):
        return f'SegmentCollection {self.segments}'

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        for segment1, segment2 in zip(self.segments, other.segments):
            if segment1 != segment2:
                return False
        return True


class GtfEntryError(Exception):
    pass


class GtfError(Exception):
    pass


class GtfEntry:
    """Represents a single GTF entry.

    Attributes:
        PARSER: Static attribute that contains a compiled regex. Used to parse a GTF line.
        ATTRIBUTE_PARSER: Static attribute that contains a compiled regex.
            Used to parse GTF entry attributes.
        _line: Raw GTF line; for internal use only. Use :attr:`line` instead.
        _chromosome: Chromosome; for internal use only. Use :attr:`chromosome` instead.
        _feature: Feature; for internal use only. Use :attr:`feature` instead.
        _start: Start; for internal use only. Use :attr:`start` instead.
        _end: End; for internal use only. Use :attr:`end` instead.
        _strand: Strand; for internal use only. Use :attr:`strand` instead.
        _attribute_str: Raw GTF entry attribute string; for internal use only.
            Use :attr:`attributes` instead.
    """
    PARSER = re.compile(
        r'''
        ^(?P<chromosome>.+?)\s+ # chromosome
        .*?\t                   # source
        (?P<feature>.+?)\s+     # feature: transcript, exon, etc.
        (?P<start>[0-9]+?)\s+   # start position (1-indexed)
        (?P<end>[0-9]+?)\s+     # end position (1-indexed, inclusive)
        .*?\s+                  # score
        (?P<strand>\+|-|\.)\s+  # +, -, . indicating strand
        .*?\s+                  # frame
        (?P<attributes>.*)      # attributes
    ''', re.VERBOSE
    )
    ATTRIBUTE_PARSER = re.compile(r'(?P<key>\S+?)\s*"(?P<value>.+?)";?')

    def __init__(self, line: str):
        """
        Args:
            line: Raw GTF line.
        """
        match = self.PARSER.match(line)
        if not match:
            raise GtfEntryError(f'Failed to parse GTF line {line}')

        self._line = line
        self._chromosome = match['chromosome']
        self._feature = match['feature']
        self._start = int(match['start'])  # 1-indexed
        self._end = int(match['end'])  # inclusive
        self._strand = match['strand']
        self._attribute_str = match['attributes']

    @property
    def line(self) -> str:
        """Raw GTF line"""
        return self._line

    @property
    def chromosome(self) -> str:
        """Chromosome"""
        return self._chromosome

    @property
    def feature(self) -> str:
        """Feature"""
        return self._feature

    @property
    def start(self) -> int:
        """Start, 1-indexed"""
        return self._start

    @property
    def end(self) -> int:
        """End, 1-indexed, inclusive"""
        return self._end

    @property
    def strand(self) -> str:
        """Strand"""
        return self._strand

    @property
    def attributes(self) -> Dict[str, str]:
        """Dictionary of attributes"""
        return dict(
            self.ATTRIBUTE_PARSER.findall(self._attribute_str.replace(' ', ''))
        )

    def to_segment(self) -> Segment:
        """Convert this GTF entry into a :class:`Segment`.

        Returns:
            The new segment
        """
        return Segment(self.start - 1, self.end)


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


def parse_gtf(
    gtf_path: str,
    features: list = ['exon', 'transcript', 'gene']
) -> Generator[GtfEntry, None, None]:
    """Parse GTF and yield only the specified features as :class:`GtfEntry` instances.

    Args:
        gtf_path: path to GTF file
        features: list of GTF features to extract

    Yields:
        GTF entries
    """
    with Gtf(gtf_path, 'r') as f:
        for entry in f:
            if entry.feature in features:
                yield entry


def genes_and_transcripts_from_gtf(gtf_path: str,
                                   use_version: bool = False
                                   ) -> Tuple[dict, dict]:
    """Parse GTF for gene and transcript information. Also, compute the introns of
    each transcript.

    Args:
        gtf_path: path to GTF file
        use_version: whether or not to use gene and transcript versions

    Returns:
        Dictionary containing gene information
        Dictionary containing transcript information
    """
    gene_infos = {}
    transcript_exons = {}
    transcript_infos = {}
    renamed = set()

    for entry in parse_gtf(gtf_path, ['exon', 'transcript', 'gene']):
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
                    '`{gene_id}-{transcript_id}`.'
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

        introns = SegmentCollection()
        exons = transcript_exons.get(transcript_id, SegmentCollection())
        if exons:
            if exons[0].start > transcript_interval.start:
                introns.add_segment(
                    Segment(transcript_interval.start, exons[0].start)
                )

            for i in range(len(exons) - 1):
                introns.add_segment(Segment(exons[i].end, exons[i + 1].start))

            if exons[-1].end < transcript_interval.end:
                introns.add_segment(
                    Segment(exons[-1].end, transcript_interval.end)
                )
        else:
            logger.warning(
                f'Gene `{attributes["gene_id"]}` transcript `{transcript_id}` has no exons. '
                'The entire transcript will be marked as an intron.'
            )
            introns.add_segment(transcript_interval)

        attributes['exons'] = exons
        attributes['introns'] = introns

    logger.debug(
        f'Parsed {len(gene_infos)} genes and {len(transcript_infos)} transcripts'
    )
    return gene_infos, transcript_infos
