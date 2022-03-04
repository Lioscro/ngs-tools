import bisect
import copy
from itertools import product
from typing import List, Optional, Set, Union

from .Segment import Segment


class SegmentCollectionError(Exception):
    pass


class SegmentCollection:
    """Class to represent a collection of integer interval segments, zero-indexed.
    The segments are always sorted and overlaps are collapsed.

    Attributes:
        _segments: List of :class:`Segment` instances; for internal use only.
            Use :attr:`segments` instead.
    """

    def __init__(self, segments: Optional[List[Segment]] = None):
        """
        Args:
            segments: List of segments. Defaults to None.
        """
        self._segments = sorted(segments) if segments else []
        if self._segments:
            self.collapse()

    @property
    def segments(self) -> List[Segment]:
        """Get a list of :class:`Segment` instances"""
        return copy.deepcopy(self._segments)

    @property
    def start(self) -> int:
        """Leftmost value of all segments. 0 if there are no segments."""
        return self._segments[0].start if self._segments else 0

    @property
    def end(self) -> int:
        """Rightmost value of all segments. 0 if there are no segments."""
        return self._segments[-1].end if self._segments else 0

    def add_segment(self, segment: Segment):
        """Add a segment to the collection.

        Args:
            segment: Segment to add
        """
        bisect.insort_left(self._segments, segment)
        self.collapse()

    def add_collection(self, collection: 'SegmentCollection'):
        """Add a collection of segments to the collection.

        Args:
            collection: Collection to add
        """
        self._segments = sorted(self._segments + collection.segments)
        self.collapse()

    def __iter__(self):
        return iter(self.segments)

    def __getitem__(self, i):
        return self._segments[i]

    def __len__(self):
        return len(self._segments)

    def __bool__(self):
        return bool(self._segments)

    def __eq__(self, other: 'SegmentCollection'):
        return self._segments == other._segments

    def invert(self, bounds: Segment) -> 'SegmentCollection':
        """Invert this :class:`SegmentCollection` within the given ``bounds``.

        Args:
            bounds: The bounds to invert with respect to.

        Returns:
            A new :class:`SegmentCollection` that is inverted

        Raises:
            SegmentCollectionError: If ``bounds`` does not entirely contain this
                collection
        """
        # First, make sure the bounds entirely contains this collection.
        if not bounds.is_superset(Segment(self.start, self.end)):
            raise SegmentCollectionError(
                'Bounds must entirely contain this collection'
            )
        # If this collection is empty, the inverse is just the bounds.
        if not self._segments:
            return SegmentCollection([Segment(bounds.start, bounds.end)])

        segments = []
        # Deal with start
        if self._segments[0].start > bounds.start:
            segments.append(Segment(bounds.start, self._segments[0].start))

        for i in range(len(self._segments) - 1):
            segments.append(
                Segment(self._segments[i].end, self._segments[i + 1].start)
            )

        # Deal with end
        if self._segments[-1].end < bounds.end:
            segments.append(Segment(self._segments[-1].end, bounds.end))

        return SegmentCollection(segments)

    def collapse(self):
        """Collapse the segments in this collection such that there are no overlapping
        segments. Any overlapping segments are merged into a single large segment.
        """
        segments = []
        combined = None
        for segment in self._segments:
            # Ignore any length 0 segments.
            if segment.width == 0:
                continue
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

        self._segments = segments

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
        for cls_segment in self._segments:
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
        segments = [segment.flank(l, left, right) for segment in self._segments]
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
        return f'{self.__class__.__name__} {self._segments}'
