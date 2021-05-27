import os
from unittest import TestCase

from joblib import delayed

from ngs_tools import utils

from . import mixins


class TestUtils(mixins.TestMixin, TestCase):

    def test_ParallelWithProgress(self):
        utils.ParallelWithProgress(
            delayed(mixins.dummy_function)() for _ in range(10)
        )

    def test_is_gzip(self):
        self.assertTrue(utils.is_gzip(self.fastq_gz_path))
        self.assertFalse(utils.is_gzip(self.fastq_path))

    def test_mkstemp(self):
        path = utils.mkstemp()
        self.assertTrue(os.path.exists(path))

    def test_mkstemp_delete(self):
        path = utils.mkstemp(delete=True)
        self.assertFalse(os.path.exists(path))
