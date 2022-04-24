import gzip
import os
import shutil
import tempfile
from unittest import TestCase, mock


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


def files_equal(file1, file2, gzipped=False):
    open_f = gzip.open if gzipped else open
    with open_f(file1, 'r') as f1, open_f(file2, 'r') as f2:
        return f1.read() == f2.read()


class TestMixin(TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @classmethod
    def setUpClass(cls):
        cls.base_dir = os.path.dirname(os.path.abspath(__file__))
        cls.fixtures_dir = os.path.join(cls.base_dir, 'fixtures')

        cls.fasta_dir = os.path.join(cls.fixtures_dir, 'fasta')
        cls.fasta_path = os.path.join(cls.fasta_dir, 'small.fa')
        cls.fasta2_path = os.path.join(cls.fasta_dir, 'not_sorted.fa')
        cls.cdna_fasta_path = os.path.join(cls.fasta_dir, 'cdna.fa')
        cls.intron_fasta_path = os.path.join(cls.fasta_dir, 'intron.fa')

        cls.fastq_dir = os.path.join(cls.fixtures_dir, 'fastq')
        cls.fastq_path = os.path.join(cls.fastq_dir, 'small.fastq')
        cls.fastq2_path = os.path.join(cls.fastq_dir, 'small2.fastq')
        cls.fastq_gz_path = os.path.join(cls.fastq_dir, 'small.fastq.gz')
        cls.fastq_paths = [
            os.path.join(cls.fastq_dir, '10xv3_1.fastq.gz'),
            os.path.join(cls.fastq_dir, '10xv3_2.fastq.gz'),
        ]
        cls.fastq2_paths = [
            os.path.join(cls.fastq_dir, 'slideseq2_1.fastq.gz'),
            os.path.join(cls.fastq_dir, 'slideseq2_2.fastq.gz'),
        ]

        cls.bam_dir = os.path.join(cls.fixtures_dir, 'bam')
        cls.bam_path = os.path.join(cls.bam_dir, 'small.bam')
        cls.bam2_path = os.path.join(cls.bam_dir, 'small2.bam')
        cls.paired_bam_path = os.path.join(cls.bam_dir, 'paired.bam')

        cls.sequence_dir = os.path.join(cls.fixtures_dir, 'sequence')
        cls.sequences_path = os.path.join(cls.sequence_dir, 'sequences.txt')
        cls.qualities_path = os.path.join(cls.sequence_dir, 'qualities.txt')

        cls.gtf_dir = os.path.join(cls.fixtures_dir, 'gtf')
        cls.gtf_path = os.path.join(cls.gtf_dir, 'not_sorted.gtf')
        cls.zero_length_gtf_path = os.path.join(cls.gtf_dir, 'zero_length.gtf')
