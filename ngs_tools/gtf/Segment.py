from typing import Optional


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
            SegmentError: If ``end < start`` or ``start < 0``
        """
        if end < start or start < 0:
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
        return f'{self.__class__.__name__} {(self.start, self.end)}'
