import copy
import itertools
import os
from abc import ABC
from enum import Enum
from typing import Dict, List, Optional, Tuple, Union

from ..fastq.Read import Read

WHITELISTS_DIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), 'whitelists'
)


class SubSequenceDefinitionError(Exception):
    pass


class SubSequenceDefinition:
    """Definition of a subsequence. This class is used to parse a subsequence out from
    a list of sequences.

    TODO: anchoring

    Attributes:
        _index: Sequence index to use (from a list of sequences); for internal use only.
            Use :attr:`index` instead.
        _start: Starting position of the subsequence; for internal use only.
            Use :attr:`start` instead.
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
            SubSequenceDefinitionError: if only one of ``start`` or ``length`` is provided
        """
        # If length is provided, start must be provided.
        if length is not None and start is None:
            raise SubSequenceDefinitionError(
                '`start` must be provided if `length` is provided'
            )
        if length is not None and length < 1:
            raise SubSequenceDefinitionError('`length` must be greater than 0')

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

    def is_overlapping(self, other: 'SubSequenceDefinition') -> bool:
        """Whether this subsequence overlaps with another subsequence.

        Args:
            other: The other :class:`SubSequenceDefinition` instance to compare to

        Returns:
            True if they are overlapping. False otherwise.
        """
        if self.index != other.index:
            return False

        # Same index
        if self.start is None or other.start is None:  # This is the whole string
            return True

        # Neither of the start positions are None
        if self.end is None and other.end is None:
            return True
        elif self.end is not None and other.end is None:
            return self.end > other.start
        elif self.end is None and other.end is not None:
            return self.start < other.end
        else:  # Neither ends are None
            return not (self.start >= other.end or self.end <= other.start)

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

    def __eq__(self, other: 'SubSequenceDefinition'):
        return (self.index, self.start,
                self.end) == (other.index, other.start, other.end)

    def __repr__(self):
        return f'{self.__class__.__name__} {(self.index, self.start, self.end)}'

    def __str__(self):
        return f'{self.index},{self.start},{self.end}'


class SubSequenceParserError(Exception):
    pass


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

    @property
    def definitions(self):
        return copy.deepcopy(self._definitions)

    def is_overlapping(self, other: 'SubSequenceParser') -> bool:
        """Whether this parser overlaps with another parser. Checks all pairwise
        combinations and returns True if any two :class:`SubSequenceDefinition`
        instances overlap.

        Args:
            other: The other :class:`SubSequenceParser` instance to compare to

        Returns:
            True if they are overlapping. False otherwise.
        """
        for parser1, parser2 in itertools.product(self._definitions,
                                                  other._definitions):
            if parser1.is_overlapping(parser2):
                return True
        return False

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
        reads: List[Read],
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

    def __eq__(self, other: 'SubSequenceParser'):
        """Check whether this parser equals another. The order of definitions
        must also be equal.
        """
        return (
            len(self._definitions) == len(other._definitions) and all(
                def1 == def2
                for def1, def2 in zip(self._definitions, other._definitions)
            )
        )

    def __iter__(self):
        return iter(self.definitions)

    def __len__(self):
        return len(self._definitions)

    def __getitem__(self, i):
        return copy.deepcopy(self._definitions[i])

    def __repr__(self):
        return f'{self.__class__.__name__} {self._definitions}'

    def __str__(self):
        return ' '.join(str(definition) for definition in self._definitions)


class ChemistryError(Exception):
    pass


class Chemistry(ABC):
    """Base class to represent any kind of chemistry.

    Attributes:
        _name: Chemistry name; for internal use only. Use :attr:`name` instead.
        _description: Chemistry description; for internal use only. Use
            :attr:`description` instead.
        _files: Dictionary containing files related to this chemistry. For internal
            use only.
    """

    def __init__(
        self,
        name: str,
        description: str,
        files: Optional[Dict[str, str]] = None
    ):
        """
        Args:
            name: Chemistry name
            description: Chemistry description
            files: Dictionary of files

        Raises:
            ChemistryError: If any of the provided files do not exist
        """
        self._name = name
        self._description = description
        self._files = files or {}

        # Check that all files exist
        for path in self._files.values():
            if not os.path.isfile(path):
                raise ChemistryError(f'File {path} does not exist')

    @property
    def name(self) -> str:
        """Chemistry name"""
        return self._name

    @property
    def description(self) -> str:
        """Chemistry description"""
        return self._description

    def has_file(self, name: str) -> bool:
        """Whether :attr:`_files` contains a file with the specified name"""
        return name in self._files

    def get_file(self, name: str) -> bool:
        """Get a file path by its name"""
        return self._files[name]


class SequencingStrand(Enum):
    UNSTRANDED = 0
    FORWARD = 1
    REVERSE = -1


class SequencingChemistry(Chemistry):
    """Base class to represent a sequencing chemistry.
    """

    def __init__(
        self, n: int, strand: SequencingStrand,
        parsers: Dict[str, SubSequenceParser], *args, **kwargs
    ):
        super().__init__(*args, **kwargs)
        self._n = n
        self._parsers = parsers
        self._strand = strand

    @property
    def n(self) -> int:
        """Number of sequences to parse at once"""
        return self._n

    @property
    def parsers(self) -> Dict[str, SubSequenceParser]:
        """Retrieve a copy of the :attr:`_parsers` dictionary."""
        return copy.deepcopy(self._parsers)

    @property
    def strand(self) -> SequencingStrand:
        """Retrieve the strandedness of the chemistry."""
        return self._strand

    @property
    def lengths(self) -> Tuple[int, ...]:
        """The expected length for each sequence, based on :attr:`parsers`.
        `None` indicates any length is expected.
        """
        l = [None] * self.n  # noqa: E741
        for parser in self._parsers.values():
            for definition in parser.definitions:
                i = definition.index
                _l = definition.end
                if _l is not None:
                    if l[i] is None:
                        l[i] = _l
                    else:
                        l[i] = max(l[i], _l)
        return tuple(l)

    def get_parser(self, name: str) -> SubSequenceParser:
        """Get a :class:`SubSequenceParser` by its name"""
        return self._parsers[name]

    def has_parser(self, name: str) -> bool:
        """Whether :attr:`_parsers` contains a parser with the specified name"""
        return name in self._parsers

    def reorder(self, reordering: List[int]) -> 'Chemistry':
        """Reorder the file indices according to the ``reordering`` list. This
        list reorders the file at each index to the value at that index.

        Args:
            reordering: List containing how to reorder file indices, where the
                file at index ``i`` of this index will now be at index
                ``reordering[i]``.

        Returns:
            A new :class:`Chemistry` instance (or the subclass)
        """
        reordered = copy.deepcopy(self)

        for parser in reordered._parsers.values():
            for _def in parser._definitions:
                _def._index = reordering[_def._index]
        reordered._n = max(reordering) + 1
        return reordered

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
        reads: List[Read],
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

    def __eq__(self, other: 'Chemistry'):
        """Check the equality of two chemistries by comparing each parser."""
        return self.n == other.n and self._parsers == other._parsers

    def __str__(self):
        return self.name

    def __repr__(self):
        return f'{self.__class__.__name__} {self.name} {self._parsers}'

    @property
    def has_barcode(self) -> bool:
        return False

    @property
    def barcode_parser(self) -> SubSequenceParser:
        raise NotImplementedError()

    @property
    def has_umi(self) -> bool:
        return False

    @property
    def umi_parser(self) -> SubSequenceParser:
        raise NotImplementedError()

    @property
    def has_whitelist(self) -> bool:
        return False

    @property
    def whitelist_path(self):
        raise NotImplementedError()

    def to_kallisto_bus_arguments(self) -> Dict[str, str]:
        """Convert this spatial chemistry definition to arguments that
        can be used as input to kallisto bus. https://www.kallistobus.tools/

        Returns:
            A Dictionary of arguments-to-value mappings. For this particular
            function, the dictionary has a single `-x` key and the value is
            a custom technology definition string, as specified in the
            kallisto manual.
        """
        if not self.has_barcode or not self.has_umi:
            raise ChemistryError(
                'Kallisto bus arguments require both `spot_barcode` and `umi` to be present.'
            )

        barcodes = []
        for _def in self.barcode_parser:
            index = _def.index
            start = _def.start or 0
            end = _def.end or 0
            barcodes.append(f'{index},{start},{end}')

        umis = []
        for _def in self.umi_parser:
            index = _def.index
            start = _def.start or 0
            end = _def.end or 0
            umis.append(f'{index},{start},{end}')

        cdnas = []
        for _def in self.cdna_parser:
            index = _def.index
            start = _def.start or 0
            end = _def.end or 0
            cdnas.append(f'{index},{start},{end}')

        return {
            '-x': f'{",".join(barcodes)}:{",".join(umis)}:{",".join(cdnas)}'
        }

    def to_starsolo_arguments(self) -> Dict[str, str]:
        """Converts this spatial chemistry definition to arguments that can
        be used as input to STARsolo.
        https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

        Returns:
            A Dictionary of arguments-to-value mappings.
        """
        args = {}
        if not self.has_barcode and not self.has_umi:
            # This must be smartseq. All cDNA definitions must be the entire
            # read, and there can be at most two.
            args['--soloType'] = 'SmartSeq'
            return args

        # Otherwise, spot barcode and UMI must exist and there must be a single
        # cDNA definition that uses the entire read.
        if not self.has_barcode or not self.has_umi:
            raise ChemistryError(
                'STARsolo requires `spot_barcode` and `umi` parsers.'
            )
        # Also, barcode and UMIs must come from the same read.
        if any(self.barcode_parser[0].index != _def.index
               for _def in list(self.barcode_parser) + list(self.umi_parser)):
            raise ChemistryError(
                'STARsolo requires spot barcode and UMI to come from the same read pair.'
            )
        # Start and end positions of spot barcode and UMI must be specified.
        if any(_def.end is None
               for _def in list(self.barcode_parser) + list(self.umi_parser)):
            raise ChemistryError(
                'STARsolo requires defined lengths for spot barcode and UMI positions.'
            )

        # Determine if CB_UMI_Simple or CB_UMI_Complex. If either barcode or
        # umi has multiple definitions, we co with complex.
        # NOTE: starsolo uses 1-indexing for start positions when CB_UMI_Simple
        # but 0-indexing for CB_UMI_Complex, while end position is inclusive
        if len(self.barcode_parser) == 1 and len(self.umi_parser) == 1:
            barcode_definition = self.barcode_parser[0]
            umi_definition = self.umi_parser[0]
            args['--soloType'] = 'CB_UMI_Simple'
            args['--soloCBstart'] = barcode_definition.start + 1
            args['--soloCBlen'] = barcode_definition.length
            args['--soloUMIstart'] = umi_definition.start + 1
            args['--soloUMIlen'] = umi_definition.length
        else:
            args['--soloType'] = 'CB_UMI_Complex'
            # No anchoring is supported yet. TODO: anchoring
            args['--soloCBposition'] = [
                f'0_{_def.start}_0_{_def.end-1}' for _def in self.barcode_parser
            ]
            args['--soloUMIposition'] = [
                f'0_{_def.start}_0_{_def.end-1}' for _def in self.umi_parser
            ]

        # Add whitelist
        args['--soloCBwhitelist'
             ] = self.whitelist_path if self.has_whitelist else 'None'
        return args
