from unittest import TestCase

from ngs_tools import chemistry, fastq

from tests.mixins import TestMixin


class TestSubSequenceDefinition(TestMixin, TestCase):

    def test_init(self):
        chemistry.SubSequenceDefinition(0, 0, 1)
        chemistry.SubSequenceDefinition(0, 1, None)
        chemistry.SubSequenceDefinition(0, None, None)
        with self.assertRaises(chemistry.ChemistryError):
            chemistry.SubSequenceDefinition(0, None, 1)
            chemistry.SubSequenceDefinition(0, 0, 0)

    def test_parse(self):
        def1 = chemistry.SubSequenceDefinition(0, 0, 1)
        def2 = chemistry.SubSequenceDefinition(0, 1, None)
        def3 = chemistry.SubSequenceDefinition(0, None, None)

        self.assertEqual('a', def1.parse(['abcd', 'efgh']))
        self.assertEqual('bcd', def2.parse(['abcd', 'efgh']))
        self.assertEqual('abcd', def3.parse(['abcd', 'efgh']))


class TestSubSequenceParser(TestMixin, TestCase):

    def test_parse(self):
        def1 = chemistry.SubSequenceDefinition(0, 0, 1)
        def2 = chemistry.SubSequenceDefinition(1, 1, None)
        parser = chemistry.SubSequenceParser(def1, def2)
        self.assertEqual(('a', 'fgh'), parser.parse(['abcd', 'efgh']))
        self.assertEqual(('afgh'),
                         parser.parse(['abcd', 'efgh'], concatenate=True))

    def test_parse_reads(self):
        def1 = chemistry.SubSequenceDefinition(0, 0, 1)
        def2 = chemistry.SubSequenceDefinition(1, 1, None)
        parser = chemistry.SubSequenceParser(def1, def2)
        read1 = fastq.Read('@1', 'ACGT', 'ABCD')
        read2 = fastq.Read('@1', 'CGTA', '1234')
        self.assertEqual((('A', 'GTA'), ('A', '234')),
                         parser.parse_reads([read1, read2]))


class TestChemistry(TestMixin, TestCase):

    def test_parse(self):
        def1 = chemistry.SubSequenceDefinition(0, 0, 1)
        def2 = chemistry.SubSequenceDefinition(1, 1, None)
        parser = chemistry.SubSequenceParser(def1, def2)
        chem = chemistry.Chemistry(
            'test', 'description', 2, {'testing': parser}
        )
        self.assertEqual({'testing': ('a', 'fgh')}, chem.parse(['abcd',
                                                                'efgh']))

    def test_parse_reads(self):
        def1 = chemistry.SubSequenceDefinition(0, 0, 1)
        def2 = chemistry.SubSequenceDefinition(1, 1, None)
        parser = chemistry.SubSequenceParser(def1, def2)
        chem = chemistry.Chemistry(
            'test', 'description', 2, {'testing': parser}
        )
        read1 = fastq.Read('@1', 'ACGT', 'ABCD')
        read2 = fastq.Read('@1', 'CGTA', '1234')
        self.assertEqual({'testing': (('A', 'GTA'), ('A', '234'))},
                         chem.parse_reads([read1, read2]))

    def test_10xv3(self):
        seq1 = 'abcdefghijklmnopqrstuvwxyzAB'
        seq2 = 'CDEFGHIJKLMNOPQRSTUVWXYZ'
        self.assertEqual({
            'cell_barcode': ('abcdefghijklmnop',),
            'umi': ('qrstuvwxyzAB',),
            'cdna': ('CDEFGHIJKLMNOPQRSTUVWXYZ',)
        }, chemistry._10X_V3.parse([seq1, seq2]))

    def test_slideseqv2(self):
        seq1 = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNO'
        seq2 = 'QRSTUVWXYZ'
        self.assertEqual({
            'spot_barcode': ('abcdefgh', 'ABCDEF'),
            'umi': ('GHIJKLMNO',),
            'cdna': ('QRSTUVWXYZ',)
        }, chemistry._SLIDESEQ_V2.parse([seq1, seq2]))

    def test_get_chemistry(self):
        self.assertEqual(chemistry._10X_V3, chemistry.get_chemistry('10xv3'))
        self.assertEqual(chemistry._10X_V3, chemistry.get_chemistry('10x-v3'))
        self.assertEqual(chemistry._10X_V3, chemistry.get_chemistry('10XV3'))
