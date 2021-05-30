from unittest import TestCase

from ngs_tools import chemistry, fastq

from tests.mixins import TestMixin


class TestSubSequenceDefinition(TestMixin, TestCase):

    def test_init(self):
        chemistry.SubSequenceDefinition(0, 0, 1)
        chemistry.SubSequenceDefinition(0, 1, None)
        chemistry.SubSequenceDefinition(0, None, None)
        with self.assertRaises(chemistry.SubSequenceDefinitionError):
            chemistry.SubSequenceDefinition(0, None, 1)
            chemistry.SubSequenceDefinition(0, 0, 0)

    def test_parse(self):
        def1 = chemistry.SubSequenceDefinition(0, 0, 1)
        def2 = chemistry.SubSequenceDefinition(0, 1, None)
        def3 = chemistry.SubSequenceDefinition(0, None, None)

        self.assertEqual('a', def1.parse(['abcd', 'efgh']))
        self.assertEqual('bcd', def2.parse(['abcd', 'efgh']))
        self.assertEqual('abcd', def3.parse(['abcd', 'efgh']))

    def test_eq(self):
        def1 = chemistry.SubSequenceDefinition(0, 0, 1)
        def2 = chemistry.SubSequenceDefinition(0, 1, None)
        def3 = chemistry.SubSequenceDefinition(0, 0, 1)
        self.assertFalse(def1 == def2)
        self.assertTrue(def1 == def3)


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

    def test_eq(self):
        def1 = chemistry.SubSequenceDefinition(0, 0, 1)
        def2 = chemistry.SubSequenceDefinition(1, 1, None)
        def3 = chemistry.SubSequenceDefinition(0, None, None)
        parser1 = chemistry.SubSequenceParser(def1, def2)
        parser2 = chemistry.SubSequenceParser(def2, def3)
        parser3 = chemistry.SubSequenceParser(def1, def2)
        self.assertFalse(parser1 == parser2)
        self.assertTrue(parser1 == parser3)


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

    def test_eq(self):
        def1 = chemistry.SubSequenceDefinition(0, 0, 1)
        def2 = chemistry.SubSequenceDefinition(1, 1, None)
        def3 = chemistry.SubSequenceDefinition(0, None, None)
        parser1 = chemistry.SubSequenceParser(def1, def2)
        parser2 = chemistry.SubSequenceParser(def2, def3)
        parser3 = chemistry.SubSequenceParser(def1, def2)
        chem1 = chemistry.Chemistry(
            'test', 'description', 2, {
                'testing': parser1,
                'testing2': parser2
            }
        )
        chem2 = chemistry.Chemistry(
            'test2', 'description2', 2, {
                'testing': parser1,
                'testing2': parser2
            }
        )
        chem3 = chemistry.Chemistry(
            'test3', 'description3', 3, {
                'testing': parser1,
                'testing2': parser2
            }
        )
        chem4 = chemistry.Chemistry(
            'test4', 'description4', 2, {
                'testing': parser2,
                'testing2': parser3
            }
        )
        self.assertFalse(chem1 == chem3)
        self.assertFalse(chem1 == chem4)
        self.assertTrue(chem1 == chem2)

    def test_10xv3(self):
        seq1 = 'abcdefghijklmnopqrstuvwxyzAB'
        seq2 = 'CDEFGHIJKLMNOPQRSTUVWXYZ'
        self.assertEqual({
            'cell_barcode': ('abcdefghijklmnop',),
            'umi': ('qrstuvwxyzAB',),
            'cdna': ('CDEFGHIJKLMNOPQRSTUVWXYZ',)
        },
                         chemistry.get_chemistry('10xv3').parse([seq1, seq2]))

    def test_to_kallisto_bus_arguments(self):
        chem = chemistry.get_chemistry('10xv3')
        self.assertEqual({'-x': '0,0,16:0,16,28:1,0,0'},
                         chem.to_kallisto_bus_arguments())

        chem = chemistry.get_chemistry('indropsv3')
        self.assertEqual({'-x': '0,0,8,1,0,8:1,8,14:2,0,0'},
                         chem.to_kallisto_bus_arguments())

        chem = chemistry.get_chemistry('smartseqv2')
        with self.assertRaises(chemistry.SingleCellChemistryError):
            chem.to_kallisto_bus_arguments()

    def test_to_starsolo_arguments(self):
        chem = chemistry.get_chemistry('10xv3')
        self.assertEqual({
            '--soloType': 'CB_UMI_Simple',
            '--soloCBstart': 1,
            '--soloCBlen': 16,
            '--soloUMIstart': 17,
            '--soloUMIlen': 12,
            '--soloCBwhitelist': chem.whitelist_path
        }, chem.to_starsolo_arguments())

        chem = chemistry.get_chemistry('indropsv1')
        self.assertEqual({
            '--soloType': 'CB_UMI_Complex',
            '--soloCBposition': ['0_0_0_10', '0_30_0_37'],
            '--soloUMIposition': ['0_42_0_47'],
            '--soloCBwhitelist': 'None'
        }, chem.to_starsolo_arguments())

        chem = chemistry.get_chemistry('indropsv3')
        with self.assertRaises(chemistry.SingleCellChemistryError):
            chem.to_starsolo_arguments()

    def test_slideseqv2(self):
        seq1 = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNO'
        seq2 = 'QRSTUVWXYZ'
        self.assertEqual({
            'spot_barcode': ('abcdefgh', 'ABCDEF'),
            'umi': ('GHIJKLMNO',),
            'cdna': ('QRSTUVWXYZ',)
        },
                         chemistry.get_chemistry('slideseqv2').parse([
                             seq1, seq2
                         ]))

    def test_get_chemistry(self):
        self.assertEqual(
            chemistry.get_chemistry('10xv3'), chemistry.get_chemistry('10x-v3')
        )
        self.assertEqual(
            chemistry.get_chemistry('10xv3'), chemistry.get_chemistry('10XV3')
        )
