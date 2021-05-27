import os
from typing import Dict, List, Optional, Tuple, Union

from . import fastq


class ChemistryError(Exception):
    pass


class SubSequenceDefinition:
    """Definition of a subsequence. This class is used to parse a subsequence out from
    a list of sequences.

    Attributes:
        _index: Sequence index to use (from a list of sequences); for internal use only.
            Use :attr:`index` instead.
        _start: Starting position of the subsequence; for internal use only.
            Use :attra:`start` instead.
        _length: Length of the subsequence; for internal use only. Use ``length`` instead.

    """

    def __init__(
        self,
        index: int,
        start: Optional[int] = None,
        length: Optional[int] = None
    ):
        """
        Args:
            index: Sequence index to extract the substring from
            start: Starting position of the substring. Defaults to None. If this
                is provided, ``length`` must also be provided.
            length: Length of the substring. Defaults to None. If this is provided,
                ``start`` must also be provided.

        Raises:
            ChemistryError: if only one of ``start`` or ``length`` is provided
        """
        # If length is provided, start must be provided.
        if length is not None and start is None:
            raise ChemistryError(
                '`start` must be provided if `length` is provided'
            )
        if length is not None and length < 1:
            raise ChemistryError('`length` must be greater than 0')

        self._index = index
        self._start = start
        self._length = length

    @property
    def index(self) -> int:
        """Sequence index"""
        return self._index

    @property
    def start(self) -> Optional[int]:
        """Substring starting position"""
        return self._start

    @property
    def end(self) -> Optional[int]:
        """Substring end position. None if :attr:`start` or :attr:`length` is None."""
        return self.start + self.length if self.start is not None and self.length is not None else None

    @property
    def length(self) -> Optional[int]:
        """Substring length. None if not provided on initialization."""
        return self._length

    def parse(self, s: List[str]) -> str:
        """Parse the given list of strings according to the arguments used to
        initialize this instance. If :attr:`start` and :attr:`length` was not provided, then
        this is simply the entire string at index :attr:`index`. Otherwise, the substring
        from position :attr:`start` of length :attr:`length` is extracted from the string
        at index :attr:`index`.

        Args:
            s: List of strings to parse

        Return:
            The parsed string
        """
        if self.start is None:
            return s[self.index]
        if len(s[self.index]) <= self.start:
            raise IndexError('string index out of range')
        if self.length is None:
            return s[self.index][self.start:]
        if len(s[self.index]) < self.end:
            raise IndexError('string index out of range')
        return s[self.index][self.start:self.end]

    def __repr__(self):
        return f'{self.__class__.__name__} {(self.index, self.start, self.end)}'

    def __str__(self):
        return f'{self.index},{self.start},{self.end}'


class SubSequenceParser:
    """Class that uses a collection of :class:`SubSequenceDefinition` instances to parse
    an entire subsequence from a list of strings.

    Attributes:
        _definitions: List of :class:`SubSequenceDefinition` instances; for internal
            use only.
    """

    def __init__(self, *definitions: SubSequenceDefinition):
        """
        Args:
            *definitions: :class:`SubSequenceDefinition` instances that are used to
                iteratively parse a list of sequences.
        """
        self._definitions = definitions

    def parse(self,
              sequences: List[str],
              concatenate: bool = False) -> Union[str, Tuple[str]]:
        """Iteratively constructs a full subsequence by applying each :class:`SubSequenceDefinition`
        in :attr:`_definitions` on the list of provided sequences. If ``concatenate=False``,
        then this function returns a tuple of length equal to the number of definitions.
        Each element of the tuple is a string that was parsed by each definition.
        Otherwise, all the parsed strings are concatenated into a single string.

        Args:
            sequences: List of sequences to parse
            concatenate: Whether or not to concatenate the parsed strings.
                Defaults to False.

        Returns:
            Concatenated parsed sequence (if ``concatenate=True``). Otherwise,
            a tuple of parsed strings.
        """
        sequence = []
        for definition in self._definitions:
            sequence.append(definition.parse(sequences))
        return ''.join(sequence) if concatenate else tuple(sequence)

    def parse_reads(
        self,
        reads: List[fastq.Read],
        concatenate: bool = False
    ) -> Tuple[Union[str, Tuple[str]], Union[str, Tuple[str]]]:
        """Behaves identically to :func:`parse`, but instead on a list of
        :class:`ngs_tools.fastq.Read` instances. :func:`parse` is called on the
        read sequences and qualities separately.

        Args:
            reads: List of reads to parse
            concatenate: Whether or not to concatenate the parsed strings.
                Defaults to False.

        Returns:
            Parsed sequence from read sequences
            Parsed sequence from quality sequences
        """
        sequences = []
        qualities = []
        for read in reads:
            sequences.append(read.sequence)
            qualities.append(read.qualities.string)
        return self.parse(sequences,
                          concatenate), self.parse(qualities, concatenate)

    def __repr__(self):
        return f'{self.__class__.__name__} {self._definitions}'

    def __str__(self):
        return ':'.join(str(definition) for definition in self._definitions)


class Chemistry:
    """Base class to represent a sequencing chemistry.

    Attributes:
        _name: Chemistry name; for internal use only. Use :attr:`name` instead.
        _description: Chemistry description; for internal use only. Use
            :attr:`description` instead.
        _n: Number of sequences (i.e. reads) that make up a single entry for this
            chemistry. For example, for paired-end reads this would be 2 (for each pair);
            for internal use only. Use :attr:`n` instead.
        _parsers: Dictionary containing :class:`SubSequenceParser` instances used to
            parse each group of :attr:`n` sequences. Each key represents a unique
            subsequence, such as cell barcode, UMI, etc. For internal use only.
    """

    def __init__(
        self, name: str, description: str, n: int,
        parsers: Dict[str, SubSequenceParser]
    ):
        """
        Args:
            name: Chemistry name
            description: Chemistry description
            n: Number of sequences
            parsers: Dictionary of parsers
        """
        self._name = name
        self._description = description
        self._n = n
        self._parsers = parsers

    @property
    def name(self) -> str:
        """Chemistry name"""
        return self._name

    @property
    def description(self) -> str:
        """Chemistry description"""
        return self._description

    @property
    def n(self) -> int:
        """Number of sequences to parse at once"""
        return self._n

    def has_parser(self, name: str) -> bool:
        """Whether :attr:`_parsers` contains a parser with the specified name"""
        return name in self._parsers

    def parse(self,
              sequences: List[str],
              concatenate: bool = False) -> Dict[str, Union[str, Tuple[str]]]:
        """Parse a list of strings using the parsers in :attr:`_parsers` and return
        a dictionary with keys corresponding to those in :attr:`_parsers`.

        Args:
            sequences: List of strings
            concatenate: Whether or not to concatenate the parsed strings.
                Defaults to False.

        Returns:
            Dictionary containing parsed strings

        Raises:
            ChemistryError: If the number sequences does not equal :attr:`n`
        """
        if len(sequences) != self.n:
            raise ChemistryError(
                f'{len(sequences)} provided but expected {self.n}'
            )
        parsed = {}
        for key, parser in self._parsers.items():
            parsed[key] = parser.parse(sequences, concatenate)
        return parsed

    def parse_reads(
        self,
        reads: List[fastq.Read],
        concatenate: bool = False,
        check_name: bool = True
    ) -> Dict[str, Tuple[Union[str, Tuple[str]], Union[str, Tuple[str]]]]:
        """Behaves identically to :func:`parse` but on a list of :class:`ngs_tools.fastq.Read`
        instances. The resulting dictionary contains tuple values, where the first
        element corresponds to the parsed read sequences, while the second corresponds to
        the parsed quality strings.

        Args:
            reads: List of :class:`ngs_tools.fastq.Read` instances
            concatenate: Whether or not to concatenate the parsed strings.
                Defaults to False.
            check_name: If True, raises :class:`ChemistryError` if all the reads
                do not have the same name. Defaults to True.

        Returns:
            Dictionary containing tuples of parsed read sequences and quality strings

        Raises:
            ChemistryError: If the number sequences does not equal :attr:`n`,
                or ``check_name=True`` and not all reads have the same name.
        """
        if len(reads) != self.n:
            raise ChemistryError(f'{len(reads)} provided but expected {self.n}')
        if check_name and not all(read.name == reads[0].name for read in reads):
            raise ChemistryError(
                'All reads must have the same name when `check_name=True`'
            )

        parsed = {}
        for key, parser in self._parsers.items():
            parsed[key] = parser.parse_reads(reads, concatenate)
        return parsed

    def __str__(self):
        return self.description

    def __repr__(self):
        return f'{self.__class__.__name__} {self.name} {self.parsers}'


class SingleCellChemistry(Chemistry):
    """Extends :class:`Chemistry` to be able to handle common single-cell
    chemistries.
    """

    def __init__(
        self,
        name: str,
        description: str,
        n: int,
        cdna_parser: SubSequenceParser,
        cell_barcode_parser: Optional[SubSequenceParser] = None,
        umi_parser: Optional[SubSequenceParser] = None,
        whitelist_path: Optional[str] = None,
    ):
        parsers = {'cdna': cdna_parser}
        if cell_barcode_parser is not None:
            parsers['cell_barcode'] = cell_barcode_parser
        if umi_parser is not None:
            parsers['umi'] = umi_parser

        super(SingleCellChemistry, self).__init__(name, description, n, parsers)
        self._whitelist_path = whitelist_path

    @property
    def has_cell_barcode(self) -> bool:
        """Whether the chemistry has a cell barcode"""
        return self.has_parser('cell_barcode')

    @property
    def has_umi(self) -> bool:
        """Whether the chemistry has a UMI"""
        return self.has_parser('umi')

    @property
    def has_whitelist(self) -> bool:
        """Whether the chemistry has a fixed predefined cell barcode whitelist"""
        return self._whitelist_path is not None

    @property
    def whitelist_path(self) -> Optional[str]:
        """Path to the whitelist. None if it does not exist."""
        return self._whitelist_path


class SpatialChemistry(Chemistry):
    """Extends :class:`Chemistry` to be able to handle common spatial chemistries.
    """

    def __init__(
        self,
        name: str,
        description: str,
        n: int,
        cdna_parser: SubSequenceParser,
        spot_barcode_parser: Optional[SubSequenceParser] = None,
        umi_parser: Optional[SubSequenceParser] = None,
        whitelist_path: Optional[str] = None,
    ):
        parsers = {'cdna': cdna_parser}
        if spot_barcode_parser is not None:
            parsers['spot_barcode'] = spot_barcode_parser
        if umi_parser is not None:
            parsers['umi'] = umi_parser

        super(SpatialChemistry, self).__init__(name, description, n, parsers)
        self._whitelist_path = whitelist_path

    @property
    def has_spot_barcode(self) -> bool:
        """Whether the chemistry has a spot barcode"""
        return self.has_parser('spot_barcode')

    @property
    def has_umi(self) -> bool:
        """Whether the chemistry has a UMI"""
        return self.has_parser('umi')

    @property
    def has_whitelist(self) -> bool:
        """Whether the chemistry has a fixed predefined spot barcode whitelist"""
        return self._whitelist_path is not None

    @property
    def whitelist_path(self) -> Optional[str]:
        """Path to the whitelist. None if it does not exist."""
        return self._whitelist_path


WHITELISTS_DIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), 'whitelists'
)

# Single cell chemistry definitions
_10X_V1 = SingleCellChemistry(
    name='10xv1',
    description='10x Genomics 3\' version 1',
    n=3,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(2)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 14)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(1, 0, 10)),
    whitelist_path=os.path.join(WHITELISTS_DIR, '10x_version1.txt.gz'),
)
_10X_V2 = SingleCellChemistry(
    name='10xv2',
    description='10x Genomics 3\' version 2',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 16)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 16, 10)),
    whitelist_path=os.path.join(WHITELISTS_DIR, '10x_version2.txt.gz'),
)
_10X_V3 = SingleCellChemistry(
    name='10xv3',
    description='10x Genomics 3\' version 3',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 16)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 16, 12)),
    whitelist_path=os.path.join(WHITELISTS_DIR, '10x_version3.txt.gz'),
)
_DROPSEQ = SingleCellChemistry(
    name='Drop-seq',
    description=(
        'Droplet-based single-cell RNA-seq chemistry developed by Macosko et al. 2015'
    ),
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 12)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 12, 8)),
    whitelist_path=os.path.join(WHITELISTS_DIR, '10x_version3.txt.gz'),
)
_SMARTSEQ_V2 = SingleCellChemistry(
    name='Smart-seq2',
    description=(
        'Plate-based single-cell RNA-seq chemistry developed by Ramskold et al. 2012'
    ),
    n=2,
    cdna_parser=SubSequenceParser(
        SubSequenceDefinition(0), SubSequenceDefinition(1)
    ),
)
_SMARTSEQ_V3 = SingleCellChemistry(
    name='Smart-seq3',
    description=(
        'Plate-based single-cell RNA-seq chemistry developed by Picelli et al. 2014'
    ),
    n=2,
    cdna_parser=SubSequenceParser(
        SubSequenceDefinition(0, 11, None), SubSequenceDefinition(1)
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 11, 8)),
)
_SINGLE_CELL_CHEMISTRIES = [
    _DROPSEQ, _10X_V1, _10X_V2, _10X_V3, _SMARTSEQ_V2, _SMARTSEQ_V3
]

# Spatial chemistry definitions
_SLIDESEQ_V2 = SpatialChemistry(
    name='Slide-seqV2',
    description=(
        'Spatial transcriptomics chemistry developed by Stickels et al. 2020'
    ),
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    spot_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(0, 0, 8), SubSequenceDefinition(0, 26, 6)
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 32, 9)),
)
_SPATIAL_CHEMISTRIES = [_SLIDESEQ_V2]

_CHEMISTRIES = _SINGLE_CELL_CHEMISTRIES + _SPATIAL_CHEMISTRIES


def get_chemistry(name: str):
    for chemistry in _CHEMISTRIES:
        if chemistry.name.replace('-', '').lower() == name.replace('-',
                                                                   '').lower():
            return chemistry

    raise ChemistryError(f'Chemistry `{name}` not found')
