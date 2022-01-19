import os
from typing import Dict, Optional

from .Chemistry import (
    Chemistry,
    SubSequenceDefinition,
    SubSequenceParser,
    WHITELISTS_DIR,
)


class SpatialChemistryError(Exception):
    pass


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
    def spot_barcode_parser(self) -> SubSequenceParser:
        """Get the spot barcode parser"""
        return self.get_parser('spot_barcode')

    @property
    def barcode_parser(self) -> SubSequenceParser:
        """Get the spot barcode parser"""
        return self.spot_barcode_parser

    @property
    def umi_parser(self) -> SubSequenceParser:
        """Get the UMI parser"""
        return self.get_parser('umi')

    @property
    def cdna_parser(self) -> SubSequenceParser:
        """Get the cDNA parser"""
        return self.get_parser('cdna')

    @property
    def has_spot_barcode(self) -> bool:
        """Whether the chemistry has a spot barcode"""
        return self.has_parser('spot_barcode')

    @property
    def has_barcode(self) -> bool:
        """Whether the chemistry has a spot barcode"""
        return self.has_spot_barcode

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

    def to_kallisto_bus_arguments(self) -> Dict[str, str]:
        """Convert this spatial chemistry definition to arguments that
        can be used as input to kallisto bus. https://www.kallistobus.tools/

        Returns:
            A Dictionary of arguments-to-value mappings. For this particular
            function, the dictionary has a single `-x` key and the value is
            a custom technology definition string, as specified in the
            kallisto manual.
        """
        if not self.has_spot_barcode or not self.has_umi:
            raise SpatialChemistryError(
                'Kallisto bus arguments require both `spot_barcode` and `umi` to be present.'
            )

        spot_barcodes = []
        for _def in self.spot_barcode_parser:
            index = _def.index
            start = _def.start or 0
            end = _def.end or 0
            spot_barcodes.append(f'{index},{start},{end}')

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
            '-x':
                f'{",".join(spot_barcodes)}:{",".join(umis)}:{",".join(cdnas)}'
        }

    def to_starsolo_arguments(self) -> Dict[str, str]:
        """Converts this spatial chemistry definition to arguments that can
        be used as input to STARsolo.
        https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

        Returns:
            A Dictionary of arguments-to-value mappings.
        """
        args = {}
        if not self.has_spot_barcode and not self.has_umi:
            # This must be smartseq. All cDNA definitions must be the entire
            # read, and there can be at most two.
            args['--soloType'] = 'SmartSeq'
            return args

        # Otherwise, spot barcode and UMI must exist and there must be a single
        # cDNA definition that uses the entire read.
        if not self.has_spot_barcode or not self.has_umi:
            raise SpatialChemistryError(
                'STARsolo requires `spot_barcode` and `umi` parsers.'
            )
        # Also, barcode and UMIs must come from the same read.
        if any(self.spot_barcode_parser[0].index != _def.index
               for _def in list(self.spot_barcode_parser) +
               list(self.umi_parser)):
            raise SpatialChemistryError(
                'STARsolo requires spot barcode and UMI to come from the same read pair.'
            )
        # Start and end positions of spot barcode and UMI must be specified.
        if any(_def.end is None for _def in list(self.spot_barcode_parser) +
               list(self.umi_parser)):
            raise SpatialChemistryError(
                'STARsolo requires defined lengths for spot barcode and UMI positions.'
            )

        # Determine if CB_UMI_Simple or CB_UMI_Complex. If either barcode or
        # umi has multiple definitions, we co with complex.
        # NOTE: starsolo uses 1-indexing for start positions when CB_UMI_Simple
        # but 0-indexing for CB_UMI_Complex, while end position is inclusive
        if len(self.spot_barcode_parser) == 1 and len(self.umi_parser) == 1:
            barcode_definition = self.spot_barcode_parser[0]
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
                f'0_{_def.start}_0_{_def.end-1}'
                for _def in self.spot_barcode_parser
            ]
            args['--soloUMIposition'] = [
                f'0_{_def.start}_0_{_def.end-1}' for _def in self.umi_parser
            ]

        # Add whitelist
        args['--soloCBwhitelist'
             ] = self.whitelist_path if self.has_whitelist else 'None'
        return args


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
_VISIUM = SpatialChemistry(
    name='Visium',
    description='10x Genomics Visium',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    spot_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 16)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 16, 12)),
    whitelist_path=os.path.join(WHITELISTS_DIR, 'visium_whitelist.txt.gz')
)
SPATIAL_CHEMISTRIES = [_SLIDESEQ_V2, _VISIUM]
