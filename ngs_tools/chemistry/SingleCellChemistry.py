import os
from typing import Optional

from .Chemistry import (
    WHITELISTS_DIR,
    SequencingChemistry,
    SequencingStrand,
    SubSequenceDefinition,
    SubSequenceParser,
)


class SingleCellChemistryError(Exception):
    pass


class SingleCellChemistry(SequencingChemistry):
    """Extends :class:`SequencingChemistry` to be able to handle common single-cell
    chemistries.
    """

    def __init__(
        self,
        name: str,
        description: str,
        n: int,
        strand: SequencingStrand,
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

        files = {}
        if whitelist_path is not None:
            files['whitelist'] = whitelist_path

        super().__init__(
            name=name,
            description=description,
            n=n,
            strand=strand,
            parsers=parsers,
            files=files
        )

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


# Single cell chemistry definitions
_10X_V1 = SingleCellChemistry(
    name='10xv1',
    description='10x Genomics 3\' version 1',
    n=3,
    strand=SequencingStrand.FORWARD,
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
    strand=SequencingStrand.FORWARD,
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
    strand=SequencingStrand.FORWARD,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 16)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 16, 12)),
    whitelist_path=os.path.join(
        WHITELISTS_DIR, '10x_version3_whitelist.txt.gz'
    ),
)

_10X_V3_ULTIMA = SingleCellChemistry(
    name='10xv3_Ultima',
    description='10x Genomics 3\' version 3 sequenced with Ultima',
    n=1,
    strand=SequencingStrand.FORWARD,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(0, 62, None)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 22, 16)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 38, 12)),
    whitelist_path=os.path.join(
        WHITELISTS_DIR, '10x_version3_whitelist.txt.gz'
    ),
)

_10X_FB = SingleCellChemistry(
    name='10xFBonly',
    description='10x Genomics Feature Barcoding',
    n=2,
    strand=SequencingStrand.FORWARD,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 16)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 16, 12)),
    whitelist_path=os.path.join(WHITELISTS_DIR, '10x_fb_whitelist.txt.gz')
)
_10X_ATAC = SingleCellChemistry(
    name='10xATAC',
    description='10x Genomics ATAC-seq',
    n=3,
    strand=SequencingStrand.FORWARD,
    cdna_parser=SubSequenceParser(
        SubSequenceDefinition(0), SubSequenceDefinition(1)
    ),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(2, 0, 16)),
    umi_parser=None,
    whitelist_path=os.path.join(WHITELISTS_DIR, '10x_atac_whitelist.txt.gz'),
)
_DROPSEQ = SingleCellChemistry(
    name='Drop-seq',
    description=(
        'Droplet-based single-cell RNA-seq chemistry developed by Macosko et al. 2015'
    ),
    n=2,
    strand=SequencingStrand.UNSTRANDED,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 12)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 12, 8)),
)
_CELSEQ_V1 = SingleCellChemistry(
    name='CEL-Seq',
    description='Hashimshony et al. 2012',
    n=2,
    strand=SequencingStrand.FORWARD,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 8)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 8, 4)),
)
_CELSEQ_V2 = SingleCellChemistry(
    name='CEL-Seq2',
    description='Hashimshony et al. 2016',
    n=2,
    strand=SequencingStrand.FORWARD,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 6, 6)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 6)),
)
_INDROPS_V1 = SingleCellChemistry(
    name='inDropsv1',
    description='Zilionis et al. 2017',
    n=2,
    strand=SequencingStrand.UNSTRANDED,
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
    strand=SequencingStrand.UNSTRANDED,
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
    strand=SequencingStrand.UNSTRANDED,
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
    strand=SequencingStrand.UNSTRANDED,
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
    strand=SequencingStrand.FORWARD,
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
    strand=SequencingStrand.UNSTRANDED,
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
    strand=SequencingStrand.UNSTRANDED,
    cdna_parser=SubSequenceParser(
        SubSequenceDefinition(0, 11, None), SubSequenceDefinition(1)
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 11, 8)),
)

_STORMSEQ = SingleCellChemistry(
    name='STORM-seq',
    description=(
        'Plate-based, ribo-reduced single-cell total RNA-seq chemistry developed by Johnson and Rhodes et al. 2022'
    ),
    n=2,
    strand=SequencingStrand.REVERSE,
    cdna_parser=SubSequenceParser(
        SubSequenceDefinition(0), SubSequenceDefinition(1, 14, None)
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(1, 0, 8)),
)
_SCI_FATE = SingleCellChemistry(
    name='Sci-fate',
    description=(
        'Single-cell RNA-seq chemistry for metabolic labeling developed by Cao et al. 2020'
    ),
    n=2,
    strand=SequencingStrand.UNSTRANDED,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 8, 10)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 8)),
    whitelist_path=os.path.join(WHITELISTS_DIR, 'sci_fate_whitelist.txt.gz'),
)
_BDWTA = SingleCellChemistry(
    name='BD Rhapsody',
    description=('Well-based single-cell RNA-seq chemistry by BD Biosciences'),
    n=2,
    strand=SequencingStrand.FORWARD,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    cell_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(0, 0, 9),
        SubSequenceDefinition(0, 21, 9),
        SubSequenceDefinition(0, 43, 9),
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 52, 8)),
    whitelist_path=os.path.join(WHITELISTS_DIR, 'BDWTA_whitelist.txt.gz'),
)
_SPLITSEQ = SingleCellChemistry(
    name='SPLiT-seq',
    description='Rosenberg et al. 2018',
    n=2,
    strand=SequencingStrand.UNSTRANDED,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(0)),
    cell_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(1, 10, 8),
        SubSequenceDefinition(1, 48, 8),
        SubSequenceDefinition(1, 78, 8),
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(1, 0, 10)),
)
_PLATE_SINGLE_CELL_CHEMISTRIES = [_SMARTSEQ_V2, _SMARTSEQ_V3, _BDWTA, _STORMSEQ]
_DROPLET_SINGLE_CELL_CHEMISTRIES = [
    _DROPSEQ, _10X_V1, _10X_V2, _10X_V3, _10X_V3_ULTIMA, _10X_FB, _10X_ATAC,
    _INDROPS_V1, _INDROPS_V2, _INDROPS_V3, _SURECELL, _SCI_FATE
]
_OTHER_SINGLE_CELL_CHEMISTRIES = [_CELSEQ_V1, _CELSEQ_V2, _SCRBSEQ, _SPLITSEQ]
SINGLE_CELL_CHEMISTRIES = (
    _PLATE_SINGLE_CELL_CHEMISTRIES + _DROPLET_SINGLE_CELL_CHEMISTRIES +
    _OTHER_SINGLE_CELL_CHEMISTRIES
)
