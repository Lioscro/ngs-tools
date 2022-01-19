import os
from typing import Dict, Optional

from .Chemistry import (
    Chemistry,
    SubSequenceDefinition,
    SubSequenceParser,
    WHITELISTS_DIR,
)


class SingleCellChemistryError(Exception):
    pass


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
        feature_map_path: Optional[str] = None,
    ):
        parsers = {'cdna': cdna_parser}
        if cell_barcode_parser is not None:
            parsers['cell_barcode'] = cell_barcode_parser
        if umi_parser is not None:
            parsers['umi'] = umi_parser

        files = {}
        if whitelist_path is not None:
            files['whitelist'] = whitelist_path
        if feature_map_path is not None:
            files['feature_map'] = feature_map_path

        super(SingleCellChemistry,
              self).__init__(name, description, n, parsers, files)

    @property
    def cell_barcode_parser(self) -> SubSequenceParser:
        """Get the cell barcode parser"""
        return self.get_parser('cell_barcode')

    @property
    def barcode_parser(self) -> SubSequenceParser:
        """Get the cell barcode parser"""
        return self.cell_barcode_parser

    @property
    def umi_parser(self) -> SubSequenceParser:
        """Get the UMI parser"""
        return self.get_parser('umi')

    @property
    def cdna_parser(self) -> SubSequenceParser:
        """Get the cDNA parser"""
        return self.get_parser('cdna')

    @property
    def has_cell_barcode(self) -> bool:
        """Whether the chemistry has a cell barcode"""
        return self.has_parser('cell_barcode')

    @property
    def has_barcode(self) -> bool:
        """Whether the chemistry has a cell barcode"""
        return self.has_cell_barcode

    @property
    def has_umi(self) -> bool:
        """Whether the chemistry has a UMI"""
        return self.has_parser('umi')

    @property
    def has_whitelist(self) -> bool:
        """Whether the chemistry has a fixed predefined cell barcode whitelist"""
        return self.has_file('whitelist')

    @property
    def whitelist_path(self) -> str:
        """Path to the whitelist"""
        return self.get_file('whitelist')

    @property
    def has_feature_map(self) -> bool:
        """Whether the chemistry has a feature barcode map"""
        return self.has_file('feature_map')

    @property
    def feature_map_path(self) -> str:
        """Path to the feature map"""
        return self.get_file('feature_map')

    def to_kallisto_bus_arguments(self) -> Dict[str, str]:
        """Convert this single-cell chemistry definition to arguments that
        can be used as input to kallisto bus. https://www.kallistobus.tools/

        Returns:
            A Dictionary of arguments-to-value mappings. For this particular
            function, the dictionary has a single `-x` key and the value is
            a custom technology definition string, as specified in the
            kallisto manual.
        """
        if not self.has_cell_barcode or not self.has_umi:
            raise SingleCellChemistryError(
                'Kallisto bus arguments require both `cell_barcode` and `umi` to be present.'
            )

        cell_barcodes = []
        for _def in self.cell_barcode_parser:
            index = _def.index
            start = _def.start or 0
            end = _def.end or 0
            cell_barcodes.append(f'{index},{start},{end}')

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
                f'{",".join(cell_barcodes)}:{",".join(umis)}:{",".join(cdnas)}'
        }

    def to_starsolo_arguments(self) -> Dict[str, str]:
        """Converts this single-cell chemistry definition to arguments that can
        be used as input to STARsolo.
        https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

        Returns:
            A Dictionary of arguments-to-value mappings.
        """
        args = {}
        if not self.has_cell_barcode and not self.has_umi:
            # This must be smartseq. All cDNA definitions must be the entire
            # read, and there can be at most two.
            args['--soloType'] = 'SmartSeq'
            return args

        # Otherwise, cell barcode and UMI must exist and there must be a single
        # cDNA definition that uses the entire read.
        if not self.has_cell_barcode or not self.has_umi:
            raise SingleCellChemistryError(
                'STARsolo requires `cell_barcode` and `umi` parsers.'
            )
        # Also, barcode and UMIs must come from the same read.
        if any(self.cell_barcode_parser[0].index != _def.index
               for _def in list(self.cell_barcode_parser) +
               list(self.umi_parser)):
            raise SingleCellChemistryError(
                'STARsolo requires cell barcode and UMI to come from the same read pair.'
            )
        # Start and end positions of cell barcode and UMI must be specified.
        if any(_def.end is None for _def in list(self.cell_barcode_parser) +
               list(self.umi_parser)):
            raise SingleCellChemistryError(
                'STARsolo requires defined lengths for cell barcode and UMI positions.'
            )

        # Determine if CB_UMI_Simple or CB_UMI_Complex. If either barcode or
        # umi has multiple definitions, we co with complex.
        # NOTE: starsolo uses 1-indexing for start positions when CB_UMI_Simple
        # but 0-indexing for CB_UMI_Complex, while end position is inclusive
        if len(self.cell_barcode_parser) == 1 and len(self.umi_parser) == 1:
            barcode_definition = self.cell_barcode_parser[0]
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
                for _def in self.cell_barcode_parser
            ]
            args['--soloUMIposition'] = [
                f'0_{_def.start}_0_{_def.end-1}' for _def in self.umi_parser
            ]

        # Add whitelist
        args['--soloCBwhitelist'
             ] = self.whitelist_path if self.has_whitelist else 'None'
        return args


# Single cell chemistry definitions
_10X_V1 = SingleCellChemistry(
    name='10xv1',
    description='10x Genomics 3\' version 1',
    n=3,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(2)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 14)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(1, 0, 10)),
    whitelist_path=os.path.join(
        WHITELISTS_DIR, '10x_version1_whitelist.txt.gz'
    ),
)

_10X_V2 = SingleCellChemistry(
    name='10xv2',
    description='10x Genomics 3\' version 2',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 16)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 16, 10)),
    whitelist_path=os.path.join(
        WHITELISTS_DIR, '10x_version2_whitelist.txt.gz'
    ),
)
_10X_V3 = SingleCellChemistry(
    name='10xv3',
    description='10x Genomics 3\' version 3',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 16)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 16, 12)),
    whitelist_path=os.path.join(
        WHITELISTS_DIR, '10x_version3_whitelist.txt.gz'
    ),
    feature_map_path=os.path.join(
        WHITELISTS_DIR, '10x_version3_feature_map.txt.gz'
    )
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
)
_CELSEQ_V1 = SingleCellChemistry(
    name='CEL-Seq',
    description='Hashimshony et al. 2012',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 8)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 8, 4)),
)
_CELSEQ_V2 = SingleCellChemistry(
    name='CEL-Seq2',
    description='Hashimshony et al. 2016',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 6, 6)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 6)),
)
_INDROPS_V1 = SingleCellChemistry(
    name='inDropsv1',
    description='Zilionis et al. 2017',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(0, 0, 11), SubSequenceDefinition(0, 30, 8)
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 42, 6)),
)
_INDROPS_V2 = SingleCellChemistry(
    name='inDropsv2',
    description='Zilionis et al. 2017',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(0)),
    cell_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(1, 0, 11), SubSequenceDefinition(1, 30, 8)
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(1, 42, 6)),
)
_INDROPS_V3 = SingleCellChemistry(
    name='inDropsv3',
    description='Zilionis et al. 2017',
    n=3,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(2)),
    cell_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(0, 0, 8), SubSequenceDefinition(1, 0, 8)
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(1, 8, 6)),
    whitelist_path=os.path.join(
        WHITELISTS_DIR, 'indrops_version3_whitelist.txt.gz'
    )
)
_SCRBSEQ = SingleCellChemistry(
    name='SCRB-seq',
    description='Soumillon et al. 2014',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 6)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 6, 10)),
)
_SURECELL = SingleCellChemistry(
    name='SureCell',
    description=(
        'Illumina Bio-Rad SureCell WTA 3\' with ddSEQ Single-Cell Isolator'
    ),
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(0, 0, 6), SubSequenceDefinition(0, 21, 6),
        SubSequenceDefinition(0, 42, 6)
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 51, 8)),
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
_SCI_FATE = SingleCellChemistry(
    name='Sci-fate',
    description=(
        'Single-cell RNA-seq chemistry for metabolic labeling developed by Cao et al. 2020'
    ),
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 8, 10)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 8)),
    whitelist_path=os.path.join(WHITELISTS_DIR, 'sci_fate_whitelist.txt.gz'),
)
_BDWTA = SingleCellChemistry(
    name='BD Rhapsody',
    description=('Well-based single-cell RNA-seq chemistry by BD Biosciences'),
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(0, 0, 9),
        SubSequenceDefinition(0, 21, 9),
        SubSequenceDefinition(0, 43, 9),
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 52, 8)),
)
_SPLITSEQ = SingleCellChemistry(
    name='SPLiT-seq',
    description='Rosenberg et al. 2018',
    n=2,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(0)),
    cell_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(1, 10, 8),
        SubSequenceDefinition(1, 48, 8),
        SubSequenceDefinition(1, 78, 8),
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(1, 0, 10)),
)
_PLATE_SINGLE_CELL_CHEMISTRIES = [_SMARTSEQ_V2, _SMARTSEQ_V3, _BDWTA]
_DROPLET_SINGLE_CELL_CHEMISTRIES = [
    _DROPSEQ, _10X_V1, _10X_V2, _10X_V3, _INDROPS_V1, _INDROPS_V2, _INDROPS_V3,
    _SURECELL, _SCI_FATE
]
_OTHER_SINGLE_CELL_CHEMISTRIES = [_CELSEQ_V1, _CELSEQ_V2, _SCRBSEQ, _SPLITSEQ]
SINGLE_CELL_CHEMISTRIES = (
    _PLATE_SINGLE_CELL_CHEMISTRIES + _DROPLET_SINGLE_CELL_CHEMISTRIES +
    _OTHER_SINGLE_CELL_CHEMISTRIES
)
