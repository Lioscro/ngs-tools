import copy
import os
from typing import Dict

from .Chemistry import (
    WHITELISTS_DIR,
    Chemistry,
    SequencingStrand,
    SubSequenceDefinition,
    SubSequenceParser,
)
from .SingleCellChemistry import _10X_FB, _10X_V3, SingleCellChemistry


class MultimodalChemistryError(Exception):
    pass


class MultimodalChemistry:
    """Represents any chemistry that is a combination of multiple chemistries.
    For example, 10x Multiome. Note that this is not a subclass of
    :class:`Chemistry`.

    TODO: Add properties similar to Chemistry class.
    """

    def __init__(
        self, name: str, description: str, chemistries: Dict[str, Chemistry]
    ):
        """
        Args:
            name: Chemistry name
            description: Chemistry description
            chemistries: Dictionary of chemistries. All keys must be unique when
                lowercased.

        Raises:
            MultimodalChemistryError: If any of the keys are not unique when
                lowercased
        """
        if len(set(key.lower()
                   for key in chemistries.keys())) != len(chemistries):
            raise MultimodalChemistryError(
                "`chemistries` contains duplicate keys when lowercased"
            )

        self._name = name
        self._description = description
        self._chemistries = chemistries

    @property
    def name(self) -> str:
        """Chemistry name"""
        return self._name

    @property
    def description(self) -> str:
        """Chemistry description"""
        return self._description

    @property
    def chemistries(self) -> Dict[str, Chemistry]:
        return copy.deepcopy(self._chemistries)

    def chemistry(self, key: str) -> Chemistry:
        return copy.deepcopy(self._chemistries[key.lower()])


_10X_FEATUREBARCODE = MultimodalChemistry(
    name='10xFB',
    description='10x Genomics Feature Barcoding',
    chemistries={
        'gex': _10X_V3,
        'fb': _10X_FB,
    }
)
_10X_MULTIOME = MultimodalChemistry(
    name='10xMultiome',
    description='10x Genomics Multiome',
    chemistries={
        'gex':
            SingleCellChemistry(
                name='10xMultiome_GEX',
                description='10x Genomics Multiome (GEX)',
                n=2,
                strand=SequencingStrand.FORWARD,
                cdna_parser=SubSequenceParser(SubSequenceDefinition(1)),
                cell_barcode_parser=SubSequenceParser(
                    SubSequenceDefinition(0, 0, 16)
                ),
                umi_parser=SubSequenceParser(SubSequenceDefinition(0, 16, 12)),
                whitelist_path=os.path.join(
                    WHITELISTS_DIR, '10x_multiome_gex_whitelist.txt.gz'
                )
            ),
        'atac':
            SingleCellChemistry(
                name='10xMultiome_ATAC',
                description='10x Genomics Multiome (ATAC)',
                n=3,
                strand=SequencingStrand.FORWARD,
                cdna_parser=SubSequenceParser(
                    SubSequenceDefinition(0), SubSequenceDefinition(1)
                ),
                cell_barcode_parser=SubSequenceParser(
                    SubSequenceDefinition(2, 0, 16)
                ),
                umi_parser=None,
                whitelist_path=os.path.join(
                    WHITELISTS_DIR, '10x_multiome_atac_whitelist.txt.gz'
                ),
            )
    }
)
MULTIMODAL_CHEMISTRIES = [_10X_FEATUREBARCODE, _10X_MULTIOME]
