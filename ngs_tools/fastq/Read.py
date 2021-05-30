import pysam


class ReadError(Exception):
    pass


class Quality:
    """Represents a Phred33 quality string.

    Attributes:
        _string: Raw quality string; for internal use only. Use :attr:`string` instead.
    """

    def __init__(self, qualities: str):
        """
        Args:
            qualities: Raw quality string
        """
        self._string = qualities

    @property
    def string(self) -> str:
        """Raw quality string"""
        return self._string

    @property
    def values(self) -> str:
        """List of quality values"""
        return list(pysam.qualitystring_to_array(self._string))

    @property
    def probs(self) -> str:
        """The quality values converted to probabilities of error"""
        return [10**(-q / 10) for q in self.values]

    def __getitem__(self, sl: slice) -> 'Quality':
        """Return a slice of the :class:`Quality`"""
        return Quality(self.string[sl])


class Read:
    """Class that represents a FASTQ read. Once the class is initialized, with the
    header, sequence, and quality strings, they should not be changed (hence using @property).
    Only Phred33 quality scores are supported, but there is no check for this.

    Attributes:
        _header: Raw header string, including the ``@``; for internal use only.
            Use :attr:`header` instead.
        _sequence: Raw sequence string; for internal use only. Use :attr:`sequence` instead.
        _qualities: :class:`Quality` instance representing the sequence qualities;
            for internal use only. Use :attr:`qualities` instead.

    """

    def __init__(self, header: str, sequence: str, qualities: str):
        """
        Args:
            header: Raw header string, including the ``@``
            sequence: Raw sequence string
            qualities: Raw qualities string

        Raises:
            ReadError: if the ``header`` does not start with ``@``
        """
        if not header.startswith('@'):
            raise ReadError(f'FASTQ header `{header}` does not start with `@`')

        self._header = header.strip()
        self._sequence = sequence.strip()
        self._qualities = Quality(qualities.strip())

    @property
    def header(self) -> str:
        """Raw header string"""
        return self._header

    @property
    def sequence(self) -> str:
        """Raw sequence string"""
        return self._sequence

    @property
    def qualities(self) -> str:
        """:class:`Quality` instance representing the sequence qualities"""
        return self._qualities

    @property
    def name(self) -> str:
        """Name of the sequence, which comes immediately after the ``@``"""
        return self.header[1:self.header.index(' ')
                           ] if ' ' in self.header else self.header[1:]

    @property
    def attributes(self) -> str:
        """String of read attributes. This is the substring after the first space in the header."""
        return self.header[self.header.index(' ') +
                           1:] if ' ' in self.header else ''
