import os
from unittest import mock, TestCase

import pysam

from ngs_tools import bam

from tests.mixins import TestMixin, tqdm_mock


class TestBam(TestMixin, TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestBam, cls).setUpClass()
        cls.tqdm_patch = mock.patch('ngs_tools.bam.tqdm', tqdm_mock)
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

    def test_filter_bam(self):
        path = os.path.join(self.temp_dir, 'test.bam')
        bam.filter_bam(self.bam_path, lambda al: al.query_name == 'read1', path)
        with pysam.AlignmentFile(path, 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))
        self.assertEqual(1, len(reads))
        self.assertEqual('read1', reads[0].query_name)
        self.assertEqual('ACTG', reads[0].query_sequence)
        self.assertEqual([32, 32, 32, 33], list(reads[0].query_qualities))
