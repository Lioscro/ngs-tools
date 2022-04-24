import os
from unittest import TestCase, mock

import pysam

from ngs_tools import bam
from tests.mixins import TestMixin, tqdm_mock


class TestBam(TestMixin, TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestBam, cls).setUpClass()
        cls.tqdm_patch = mock.patch('ngs_tools.bam.progress', tqdm_mock)
        cls.tqdm_patch.start()

    @classmethod
    def tearDownClass(cls):
        super(TestBam, cls).tearDownClass()
        cls.tqdm_patch.stop()

    def test_map_bam(self):
        self.assertEqual([
            'ACTG', 'AAGT'
        ], list(bam.map_bam(self.bam_path, lambda al: al.query_sequence)))

    def test_apply_bam(self):
        path = os.path.join(self.temp_dir, 'test.bam')

        def apply_func(al):
            if al.query_name == 'read1':
                al.set_tag('TS', 'a')
                return al

        bam.apply_bam(self.bam_path, apply_func, path)
        with pysam.AlignmentFile(path, 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))
        self.assertEqual(1, len(reads))
        self.assertEqual('read1', reads[0].query_name)
        self.assertEqual('ACTG', reads[0].query_sequence)
        self.assertEqual([32, 32, 32, 33], list(reads[0].query_qualities))
        self.assertEqual('a', reads[0].get_tag('TS'))

    def test_tag_bam_with_fastq(self):
        path = os.path.join(self.temp_dir, 'test.bam')

        def tag_func(read):
            return {'TS': read.sequence[:2]}

        bam.tag_bam_with_fastq(self.bam_path, self.fastq2_path, tag_func, path)
        with pysam.AlignmentFile(path, 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))

        self.assertEqual(2, len(reads))
        self.assertEqual('read1', reads[0].query_name)
        self.assertEqual('ACTG', reads[0].query_sequence)
        self.assertEqual([32, 32, 32, 33], list(reads[0].query_qualities))
        self.assertEqual('TT', reads[0].get_tag('TS'))
        self.assertEqual('read2', reads[1].query_name)
        self.assertEqual('AAGT', reads[1].query_sequence)
        self.assertEqual([14, 14, 14, 0], list(reads[1].query_qualities))
        self.assertEqual('CC', reads[1].get_tag('TS'))

    def test_tag_bam_with_multiple_fastqs(self):
        path = os.path.join(self.temp_dir, 'test.bam')

        def tag_func1(read):
            return {'T1': read.sequence[:2]}

        def tag_func2(read):
            return {'T2': read.sequence[2:]}

        bam.tag_bam_with_fastq(
            self.bam_path, [self.fastq_path, self.fastq2_path],
            [tag_func1, tag_func2], path
        )
        with pysam.AlignmentFile(path, 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))

        self.assertEqual(2, len(reads))
        self.assertEqual('read1', reads[0].query_name)
        self.assertEqual('ACTG', reads[0].query_sequence)
        self.assertEqual([32, 32, 32, 33], list(reads[0].query_qualities))
        self.assertEqual('AC', reads[0].get_tag('T1'))
        self.assertEqual('TT', reads[0].get_tag('T2'))
        self.assertEqual('read2', reads[1].query_name)
        self.assertEqual('AAGT', reads[1].query_sequence)
        self.assertEqual([14, 14, 14, 0], list(reads[1].query_qualities))
        self.assertEqual('AA', reads[1].get_tag('T1'))
        self.assertEqual('CC', reads[1].get_tag('T2'))

    def test_filter_bam(self):
        path = os.path.join(self.temp_dir, 'test.bam')
        bam.filter_bam(self.bam_path, lambda al: al.query_name == 'read1', path)
        with pysam.AlignmentFile(path, 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))
        self.assertEqual(1, len(reads))
        self.assertEqual('read1', reads[0].query_name)
        self.assertEqual('ACTG', reads[0].query_sequence)
        self.assertEqual([32, 32, 32, 33], list(reads[0].query_qualities))

    def test_split_bam_n(self):
        prefix = os.path.join(self.temp_dir, 'test')
        self.assertEqual({
            '0': (f'{prefix}_0.bam', 3),
            '1': (f'{prefix}_1.bam', 1)
        }, bam.split_bam(self.bam2_path, prefix, n=2))
        with pysam.AlignmentFile(f'{prefix}_0.bam', 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))
            self.assertEqual(3, len(reads))
            self.assertEqual(['read1', 'read2', 'read3'],
                             [read.query_name for read in reads])
        with pysam.AlignmentFile(f'{prefix}_1.bam', 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))
            self.assertEqual(1, len(reads))
            self.assertEqual(['read4'], [read.query_name for read in reads])

    def test_split_bam_func(self):
        prefix = os.path.join(self.temp_dir, 'test')
        self.assertEqual({
            'test': (f'{prefix}_test.bam', 2),
            'test2': (f'{prefix}_test2.bam', 2)
        },
                         bam.split_bam(
                             self.bam2_path,
                             prefix,
                             split_func=lambda al: al.get_tag('RG')
                         ))
        with pysam.AlignmentFile(f'{prefix}_test.bam', 'rb',
                                 check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))
            self.assertEqual(2, len(reads))
            self.assertEqual(['read1', 'read2'],
                             [read.query_name for read in reads])
        with pysam.AlignmentFile(f'{prefix}_test2.bam', 'rb',
                                 check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))
            self.assertEqual(2, len(reads))
            self.assertEqual(['read3', 'read4'],
                             [read.query_name for read in reads])

    def test_split_bam_paired(self):
        prefix = os.path.join(self.temp_dir, 'test')
        self.assertEqual({
            '0': (f'{prefix}_0.bam', 3),
            '1': (f'{prefix}_1.bam', 1)
        }, bam.split_bam(self.paired_bam_path, prefix, n=2))
        with pysam.AlignmentFile(f'{prefix}_0.bam', 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))
            self.assertEqual(3, len(reads))
            self.assertEqual(['read1', 'read2', 'read1'],
                             [read.query_name for read in reads])
        with pysam.AlignmentFile(f'{prefix}_1.bam', 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))
            self.assertEqual(1, len(reads))
            self.assertEqual(['read3'], [read.query_name for read in reads])

    def test_count_bam(self):
        self.assertEqual(2, bam.count_bam(self.bam_path))

    def test_count_bam_filter(self):
        self.assertEqual(0, bam.count_bam(self.bam_path, lambda al: False))
