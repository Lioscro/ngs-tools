from unittest import TestCase

from joblib import delayed

from ngsutils import utils

from tests.mixins import TestMixin, dummy_function


class TestUtils(TestMixin, TestCase):

    def test_ParallelWithProgress(self):
        utils.ParallelWithProgress(delayed(dummy_function)() for _ in range(10))

    def test_is_gzip(self):
        self.assertTrue(utils.is_gzip(self.fastq_gz_path))
        self.assertFalse(utils.is_gzip(self.fastq_path))
