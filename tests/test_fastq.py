import os
from unittest import mock, TestCase

import pysam

from ngs_tools import fastq

from . import mixins


class TestFastq(mixins.TestMixin, TestCase):

    def test_read(self):
        fq = fastq.Fastq(self.fastq_path, 'r')
        reads = list(fq)
        self.assertEqual(2, len(reads))

        self.assertEqual('@read1 extra tags', reads[0].header)
        self.assertEqual('read1', reads[0].name)
        self.assertEqual('extra tags', reads[0].attributes)
        self.assertEqual('ACTG', reads[0].sequence)
        self.assertEqual('AAAB', reads[0].qualities.string)
        self.assertEqual([32, 32, 32, 33], reads[0].qualities.values)

        self.assertEqual('@read2 extra tags2', reads[1].header)
        self.assertEqual('read2', reads[1].name)
        self.assertEqual('extra tags2', reads[1].attributes)
        self.assertEqual('AAGT', reads[1].sequence)
        self.assertEqual('///!', reads[1].qualities.string)
        self.assertEqual([14, 14, 14, 0], reads[1].qualities.values)

    def test_gz_read(self):
        fq = fastq.Fastq(self.fastq_gz_path, 'r')
        reads = list(fq)
        self.assertEqual(2, len(reads))

        self.assertEqual('@read1 extra tags', reads[0].header)
        self.assertEqual('read1', reads[0].name)
        self.assertEqual('extra tags', reads[0].attributes)
        self.assertEqual('ACTG', reads[0].sequence)
        self.assertEqual('AAAB', reads[0].qualities.string)
        self.assertEqual([32, 32, 32, 33], reads[0].qualities.values)

        self.assertEqual('@read2 extra tags2', reads[1].header)
        self.assertEqual('read2', reads[1].name)
        self.assertEqual('extra tags2', reads[1].attributes)
        self.assertEqual('AAGT', reads[1].sequence)
        self.assertEqual('///!', reads[1].qualities.string)
        self.assertEqual([14, 14, 14, 0], reads[1].qualities.values)

    def test_write(self):
        path = os.path.join(self.temp_dir, 'test.fastq')
        fq = fastq.Fastq(path, 'w')
        fq.write(fastq.Read('@header', 'ACTG', '////'))
        fq.close()
        with open(path, 'r') as f:
            self.assertEqual('@header\nACTG\n+\n////\n', f.read())

    def test_fastq_to_bam(self):
        path = os.path.join(self.temp_dir, 'test.bam')
        with mock.patch('ngs_tools.fastq.tqdm', mixins.tqdm_mock):
            fastq.fastq_to_bam(self.fastq_path, path, 'read_group')

        with pysam.AlignmentFile(path, 'rb', check_sq=False) as f:
            reads = list(f.fetch(until_eof=True))

        self.assertEqual(2, len(reads))
        self.assertEqual('read1', reads[0].query_name)
        self.assertEqual('ACTG', reads[0].query_sequence)
        self.assertEqual([32, 32, 32, 33], list(reads[0].query_qualities))
        self.assertEqual('read_group', reads[0].get_tag('RG'))
        self.assertEqual('read2', reads[1].query_name)
        self.assertEqual('AAGT', reads[1].query_sequence)
        self.assertEqual([14, 14, 14, 0], list(reads[1].query_qualities))
        self.assertEqual('read_group', reads[1].get_tag('RG'))
