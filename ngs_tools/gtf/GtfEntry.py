import re
from typing import Dict

from .Segment import Segment


class GtfEntryError(Exception):
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
