import re

from .Chemistry import (
    Chemistry,
    ChemistryError,
    SequencingChemistry,
    SequencingStrand,
    SubSequenceDefinition,
    SubSequenceDefinitionError,
    SubSequenceParser,
    SubSequenceParserError,
)
from .MultimodalChemistry import (
    MULTIMODAL_CHEMISTRIES,
    MultimodalChemistry,
    MultimodalChemistryError,
)
from .SingleCellChemistry import (
    SINGLE_CELL_CHEMISTRIES,
    SingleCellChemistry,
    SingleCellChemistryError,
)
from .SpatialChemistry import (
    SPATIAL_CHEMISTRIES,
    SpatialChemistry,
    SpatialResolution,
    SpatialSequencingChemistry,
)

VERSION_PARSER = re.compile(r'v?\d+$')
CHEMISTRIES = SINGLE_CELL_CHEMISTRIES + SPATIAL_CHEMISTRIES + MULTIMODAL_CHEMISTRIES


def _clean_name(name: str):
    """Internal helper function to clean chemistry names.

    Args:
        name: String name of the chemistry.

    Returns:
        Cleaned name
    """
    name = name.lower().replace('-', '').replace(' ', '')
    version = 1
    base_name = name

    version_search = VERSION_PARSER.search(name)
    if version_search:
        version_suffix = name[version_search.start(0):]
        version = int(
            version_suffix[1:] if version_suffix[0] == 'v' else version_suffix
        )

        base_name = name[:version_search.start(0)]

    # Cleaned name is of the form {base_name}-{version}.
    return f'{base_name}-{version}'


def get_chemistry(name: str):
    """Fetch a :class:`Chemistry` definition by name. Uses some regex magic to
    correctly deal with chemistry versioning at the end of the name. For instance,
    ``10x2`` is interpreted the same as ``10xv2``.

    See :mod:`.SingleCellChemistry` and :mod:`.SpatialChemistry` for available
    chemistries.

    Args:
        name: String name of the chemistry. Any dashes (`-`) or capitalization
            are ignored.

    Returns:
        The matching chemistry.

    Raises:
        ChemistryError: If the chemistry could not be found.
    """
    cleaned_name = _clean_name(name)

    for chemistry in CHEMISTRIES:
        if cleaned_name == _clean_name(chemistry.name):
            return chemistry

    raise ChemistryError(f'Chemistry `{name}` not found')
