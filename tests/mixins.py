import os
import shutil
import tempfile
from unittest import mock, TestCase


def tqdm_mock(*args, **kwargs):
    if len(args) > 0:
        iterable = args[0]
        try:
            iter(iterable)
        except TypeError:
            return mock.MagicMock()
        return iterable
    return mock.MagicMock()


def dummy_function(*args, **kwargs):
    return mock.MagicMock()


class TestMixin(TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @classmethod
    def setUpClass(cls):
        cls.base_dir = os.path.dirname(os.path.abspath(__file__))
        cls.fixtures_dir = os.path.join(cls.base_dir, 'fixtures')
        cls.fastq_dir = os.path.join(cls.fixtures_dir, 'fastq')
        cls.fastq_path = os.path.join(cls.fastq_dir, 'small.fastq')
        cls.fastq2_path = os.path.join(cls.fastq_dir, 'small2.fastq')
        cls.fastq_gz_path = os.path.join(cls.fastq_dir, 'small.fastq.gz')

        cls.bam_dir = os.path.join(cls.fixtures_dir, 'bam')
        cls.bam_path = os.path.join(cls.bam_dir, 'small.bam')

        cls.sequence_dir = os.path.join(cls.fixtures_dir, 'sequence')
        cls.sequences_path = os.path.join(cls.sequence_dir, 'sequences.txt')
        cls.qualities_path = os.path.join(cls.sequence_dir, 'qualities.txt')
