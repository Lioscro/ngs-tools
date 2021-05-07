import warnings

import pysam

from ngs_tools.logging import silence_logger

pysam.set_verbosity(0)
warnings.filterwarnings('ignore')
silence_logger('ngs_tools')
