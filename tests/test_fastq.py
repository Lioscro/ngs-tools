import os
from unittest import TestCase, mock

import pysam

from ngs_tools import chemistry, fastq

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
        with mock.patch('ngs_tools.fastq.progress', mixins.tqdm_mock):
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

    def test_fastqs_to_bam_with_chemistry_10xv3(self):
        bam_path = fastq.fastqs_to_bam_with_chemistry(
            self.fastq_paths,
            chemistry.get_chemistry('10xv3'), {
                'umi': ('UR', 'UY'),
                'cell_barcode': ('CR', 'CY')
            },
            os.path.join(self.temp_dir, 'test.bam'),
            name='test'
        )
        with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as f:
            alignments = list(f.fetch(until_eof=True))
        self.assertEqual(2, len(alignments))
        self.assertEqual([
            'M03718:773:000000000-JKHP3:1:1101:18272:1693',
            'M03718:773:000000000-JKHP3:1:1101:17963:1710',
        ], [al.query_name for al in alignments])
        self.assertEqual([
            read.sequence for read in fastq.Fastq(self.fastq_paths[1])
        ], [al.query_sequence for al in alignments])
        self.assertEqual([
            read.qualities.string for read in fastq.Fastq(self.fastq_paths[1])
        ], [
            pysam.array_to_qualitystring(al.query_qualities)
            for al in alignments
        ])
        self.assertEqual({
            ('UR', 'CCAAAACAGTTT'),
            ('UY', 'CEE0C0BA0DFG'),
            ('CR', 'TACGTCATCTCCTACG'),
            ('CY', '1111AFAFFFBFGGFE'),
            ('RG', 'test'),
        }, set(alignments[0].get_tags()))
        self.assertEqual({
            ('UR', 'ATTCCTGAGTCA'),
            ('UY', 'BFGFGFF10FG1'),
            ('CR', 'TTAGATCGTTAGAAAC'),
            ('CY', '1>>11DFAFAAAFFGG'),
            ('RG', 'test'),
        }, set(alignments[1].get_tags()))

    def test_fastqs_to_bam_with_chemistry_slideseq2(self):
        bam_path = fastq.fastqs_to_bam_with_chemistry(
            self.fastq2_paths,
            chemistry.get_chemistry('slideseq2'), {
                'umi': ('UR', 'UY'),
                'spot_barcode': ('CR', 'CY')
            },
            os.path.join(self.temp_dir, 'test.bam'),
            name='test'
        )
        with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as f:
            alignments = list(f.fetch(until_eof=True))
        self.assertEqual(2, len(alignments))
        self.assertEqual([
            'NB501583:801:H7JLTBGXH:1:11101:20912:1050',
            'NB501583:801:H7JLTBGXH:1:11101:8670:1050',
        ], [al.query_name for al in alignments])
        self.assertEqual([
            read.sequence for read in fastq.Fastq(self.fastq2_paths[1])
        ], [al.query_sequence for al in alignments])
        self.assertEqual([
            read.qualities.string for read in fastq.Fastq(self.fastq2_paths[1])
        ], [
            pysam.array_to_qualitystring(al.query_qualities)
            for al in alignments
        ])
        self.assertEqual({
            ('UR', 'TTTTTTTTT'),
            ('UY', 'EEEEEEEEE'),
            ('CR', 'CTTTGNTCAATGTT'),
            ('CY', 'AAAAA#EEAEEEEE'),
            ('RG', 'test'),
        }, set(alignments[0].get_tags()))
        self.assertEqual({
            ('UR', 'AGTGTCTCA'),
            ('UY', 'EAEAEAEEE'),
            ('CR', 'CTCTTNATCCTCAT'),
            ('CY', 'AAAAA#EEE/EAE/'),
            ('RG', 'test'),
        }, set(alignments[1].get_tags()))
