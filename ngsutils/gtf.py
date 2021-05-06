import bisect
import re
from itertools import product

from . import utils


class Segment:
    """Class to represent an integer interval segment, zero-indexed.

    :param start: start position
    :type start: int
    :param end: end position
    :type end: int
    """

    def __init__(self, start, end):
        assert end > start
        self._start = start
        self._end = end

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def width(self):
        return self.end - self.start

    def is_in(self, i):
        return i >= self.start and i < self.end

    def is_exclusive(self, segment):
        return self.end <= segment.start or self.start >= segment.end

    def is_overlapping(self, segment):
        return not self.is_exclusive(segment)

    def is_subset(self, segment):
        return self.start >= segment.start and self.end <= segment.end

    def is_superset(self, segment):
        return self.start <= segment.start and self.end >= segment.end

    def __iter__(self):
        return iter((self.start, self.end))

    def __eq__(self, other):
        return (self.start, self.end) == (other.start, other.end)

    def __lt__(self, other):
        return (self.start, self.end) < (other.start, other.end)

    def __gt__(self, other):
        return (self.start, self.end) > (other.start, other.end)

    def __str__(self):
        return str((self.start, self.end))

    def __repr__(self):
        return str(self)


class SegmentCollection:
    """Class to represent a collection of integer interval segments, zero-indexed.

    :param segments: list of initial segments, defaults to `None`
    :type segments: list, optional
    """

    def __init__(self, segments=None):
        self.segments = sorted(segments) if segments else []
        if self.segments:
            self.collapse()

    @property
    def start(self):
        return self.segments[0].start if self.segments else -1

    @property
    def end(self):
        return self.segments[-1].end if self.segments else -1

    def add_segment(self, segment):
        bisect.insort_left(self.segments, segment)
        self.collapse()

    def add_collection(self, collection):
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

    def span_is_exclusive(self, collection):
        return self.end <= collection.start or self.start >= collection.end

    def is_overlapping(self, collection):
        if self.span_is_exclusive(collection):
            return False
        return any(
            segment_1.is_overlapping(segment_2)
            for segment_1, segment_2 in product(self, collection)
        )

    def is_subset(self, collection):
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

    def is_superset(self, collection):
        return collection.is_subset(self)

    @classmethod
    def from_positions(cls, positions):
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
    def from_collections(cls, *collections):
        segments = []
        for collection in collections:
            segments.extend(collection.segments)
        return cls(segments=segments)

    def __str__(self):
        return f'SegmentCollection {str(self.segments)}'

    def __repr__(self):
        return str(self)

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
    def line(self):
        return self._line

    @property
    def chromosome(self):
        return self._chromosome

    @property
    def feature(self):
        return self._feature

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    @property
    def attributes(self):
        return dict(
            self.ATTRIBUTE_PARSER.findall(self._attribute_str.replace(' ', ''))
        )

    def to_segment(self) -> Segment:
        return Segment(self.start - 1, self.end)


class Gtf(utils.FileWrapper):

    def read(self) -> GtfEntry:
        if self.mode != 'r':
            raise GtfError(f'Can not read from file in mode `{self.mode}`')
        if self.closed:
            raise GtfError('Can not read from closed file')

        while True:
            line = next(self.fp)
            if not line.startswith('#'):
                break

        return GtfEntry(line)

    def write(self, entry: GtfEntry):
        if self.mode != 'w':
            raise GtfError(f'Can not write to file in mode `{self.mode}`')
        if self.closed:
            raise GtfError('Can not write to closed file')

        self.fp.write(f'{entry.line}\n')


# def parse_gtf(gtf_path):
#     """Parse GTF for gene and transcript informations.
#
#     :param gtf_path: path to GTF
#     :type gtf_path: str
#
#     :return: (gene information, transcript information)
#     :rtype: (dictionary, dictionary)
#     """
#     gtf = GTF(gtf_path)
#
#     gene_infos = {}
#     transcript_exons = {}
#     transcript_infos = {}
#     count = 0
#     for gtf_entry in gtf.entries():
#         count += 1
#
#         feature = gtf_entry['feature']
#         start = gtf_entry['start'] - 1
#         end = gtf_entry['end']
#         strand = gtf_entry['strand']
#         chrom = gtf_entry['seqname']
#
#         # Ignore any features that are not exon, transcript, or gene
#         if feature in ('exon', 'transcript', 'gene'):
#             # IMPORTANT: every feature must have gene_id
#             gene_id = gtf_entry['group']['gene_id']
#             gene_name = gtf_entry['group'].get('gene_name')
#
#             # Update gene info
#             gene_info = gene_infos.setdefault(
#                 gene_id, {
#                     'segment': Segment(start, end),
#                     'chr': chrom,
#                     'strand': strand,
#                     'gene_name': gene_name,
#                     'transcripts': []
#                 }
#             )
#             segment = gene_info['segment']
#             gene_info['segment'] = Segment(min(segment.start, start), max(segment.end, end))
#             gene_info['gene_name'] = gene_name or gene_info['gene_name']
#
#             # Transcript and exon features
#             # Some GTFs (specifically gencode) add transcript_id to gene features,
#             # but in these cases gene_id = transcript_id, so we ignore these explicitly.
#             if feature != 'gene':
#                 # IMPORTANT: every transcript, exon feature must have transcript_id
#                 transcript_id = gtf_entry['group'].get('transcript_id')
#                 if not transcript_id:
#                     logger.warning(f'Gene `{gene_id}` feature `{feature}` does not have a transcript_id.')
#                 if gene_id != transcript_infos.get(transcript_id, {}).get('gene_id', gene_id):
#                     logger.warning(
#                         f'Transcript `{transcript_id}` is assigned to multiple genes. '
#                         f'Renaming to `{gene_id}-{transcript_id}`.'
#                     )
#                     transcript_id = f'{gene_id}-{transcript_id}'
#
#                 gene_info['transcripts'].append(transcript_id)
#                 transcript_info = transcript_infos.setdefault(
#                     transcript_id, {
#                         'gene_id': gene_id,
#                         'chr': chrom,
#                         'segment': Segment(start, end),
#                         'strand': strand
#                     }
#                 )
#                 segment = transcript_info['segment']
#                 transcript_info['segment'] = Segment(min(segment.start, start), max(segment.end, end))
#
#                 if feature == 'exon':
#                     transcript_exons.setdefault(transcript_id, SegmentCollection()).add_segment(Segment(start, end))
#
#     # Clean gene infos so that they link to only transcripts existing in transcript_infos
#     for gene_id, attributes in gene_infos.items():
#         cleaned = list(set(t for t in attributes['transcripts'] if t in transcript_infos))
#         attributes['transcripts'] = cleaned
#         if not cleaned:
#             logger.warning(
#                 f'Gene `{gene_id}` has no transcripts. '
#                 f'The entire gene will be marked as a transcript with ID `{gene_id}` and an exon.'
#             )
#             attributes['transcripts'].append(gene_id)
#             transcript_infos[gene_id] = {
#                 'gene_id': gene_id,
#                 'chr': attributes['chr'],
#                 'segment': attributes['segment'],
#                 'strand': attributes['strand']
#             }
#             transcript_exons[gene_id] = SegmentCollection(segments=[attributes['segment']])
#
#     # Calculate introns
#     for transcript_id, attributes in transcript_infos.items():
#         transcript_interval = attributes['segment']
#
#         introns = SegmentCollection()
#         exons = transcript_exons.get(transcript_id, SegmentCollection())
#         if exons:
#             if exons[0].start > transcript_interval.start:
#                 introns.add_segment(Segment(transcript_interval.start, exons[0].start))
#
#             for i in range(len(exons) - 1):
#                 introns.add_segment(Segment(exons[i].end, exons[i + 1].start))
#
#             if exons[-1].end < transcript_interval.end:
#                 introns.add_segment(Segment(exons[-1].end, transcript_interval.end))
#         else:
#             logger.warning(
#                 f'Gene `{attributes["gene_id"]}` transcript `{transcript_id}` has no exons. '
#                 'The entire transcript will be marked as an intron.'
#             )
#             introns.add_segment(transcript_interval)
#
#         attributes['exons'] = exons
#         attributes['introns'] = introns
#
#     logger.debug(f'Parsed {len(gene_infos)} genes and {len(transcript_infos)} transcripts from {count} GTF entries')
#     return gene_infos, transcript_infos
