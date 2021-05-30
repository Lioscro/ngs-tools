from typing import Optional

from .Chemistry import Chemistry, SubSequenceDefinition, SubSequenceParser


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
SPATIAL_CHEMISTRIES = [_SLIDESEQ_V2]
