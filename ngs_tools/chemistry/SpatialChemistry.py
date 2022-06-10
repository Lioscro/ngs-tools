import os
from typing import NamedTuple, Optional

from typing_extensions import Literal

from .Chemistry import (
    WHITELISTS_DIR,
    Chemistry,
    SequencingChemistry,
    SequencingStrand,
    SubSequenceDefinition,
    SubSequenceParser,
)


class SpatialResolution(NamedTuple):
    scale: float = 1.0
    unit: Optional[Literal["nm", "um", "mm"]] = None


class SpatialChemistryError(Exception):
    pass


class SpatialChemistry(Chemistry):
    """Extends :class:`Chemistry` to be able to handle spatial chemistries.
    """

    def __init__(self, resolution: SpatialResolution, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._resolution = resolution

    @property
    def resolution(self) -> SpatialResolution:
        """Get the spatial resolution as a :class:`SpatialResolution` object.
        """
        return self._resolution


class SpatialSequencingChemistry(SpatialChemistry, SequencingChemistry):
    """Extends :class:`SequencingChemistry` to be able to handle common spatial chemistries.
    """

    def __init__(
        self,
        name: str,
        description: str,
        resolution: SpatialResolution,
        n: int,
        strand: SequencingStrand,
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

        files = {}
        if whitelist_path is not None:
            files['whitelist'] = whitelist_path

        super().__init__(
            name=name,
            description=description,
            resolution=resolution,
            n=n,
            strand=strand,
            parsers=parsers,
            files=files,
        )

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
        return self.has_file('whitelist')

    @property
    def whitelist_path(self) -> Optional[str]:
        """Path to the whitelist. None if it does not exist."""
        return self.get_file('whitelist')


# Spatial chemistry definitions
_SLIDESEQ_V2 = SpatialSequencingChemistry(
    name='Slide-seqV2',
    description=(
        'Spatial transcriptomics chemistry developed by Stickels et al. 2020'
    ),
    resolution=SpatialResolution(10., 'um'),
    n=2,
    strand=SequencingStrand.UNSTRANDED,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    spot_barcode_parser=SubSequenceParser(
        SubSequenceDefinition(0, 0, 8), SubSequenceDefinition(0, 26, 6)
    ),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 32, 9)),
)
_VISIUM = SpatialSequencingChemistry(
    name='Visium',
    description='10x Genomics Visium',
    resolution=SpatialResolution(55., 'um'),
    n=2,
    strand=SequencingStrand.FORWARD,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    spot_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 16)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 16, 12)),
    whitelist_path=os.path.join(WHITELISTS_DIR, 'visium_whitelist.txt.gz')
)
_STEREOSEQ = SpatialSequencingChemistry(
    name='Stereo-seq',
    description='BGI Stereo-seq',
    resolution=SpatialResolution(0.5, 'um'),
    n=2,
    strand=SequencingStrand.UNSTRANDED,
    cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
    spot_barcode_parser=SubSequenceParser(SubSequenceDefinition(0, 0, 25)),
    umi_parser=SubSequenceParser(SubSequenceDefinition(0, 25, 10)),
)
_COSMX = SpatialChemistry(
    name='CosMx',
    description='NanoString CosMx Spatial Molecular Imager',
    resolution=SpatialResolution(0.18, 'um'),
)
_SEQUENCING_SPATIAL_CHEMISTRIES = [_SLIDESEQ_V2, _VISIUM, _STEREOSEQ]
_INSITU_SPATIAL_CHEMISTRIES = [_COSMX]
SPATIAL_CHEMISTRIES = _SEQUENCING_SPATIAL_CHEMISTRIES + _INSITU_SPATIAL_CHEMISTRIES
