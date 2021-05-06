import os
from unittest import TestCase

from ngsutils import fasta

from . import mixins


class TestFasta(mixins.TestMixin, TestCase):

    def test_read(self):
        fa = fasta.Fasta(self.fasta_path, 'r')
        entries = list(fa)
        self.assertEqual(2, len(entries))

        self.assertEqual('>seq1 attr1:val1 attr2:val2', entries[0].header)
        self.assertEqual('seq1', entries[0].name)
        self.assertEqual({
            'attr1': 'val1',
            'attr2': 'val2'
        }, entries[0].attributes)
        self.assertEqual('ACTG', entries[0].sequence)

        self.assertEqual('>seq2 attr3:val3', entries[1].header)
        self.assertEqual('seq2', entries[1].name)
        self.assertEqual({'attr3': 'val3'}, entries[1].attributes)
        self.assertEqual('AACC', entries[1].sequence)

    def test_write(self):
        path = os.path.join(self.temp_dir, 'test.fa')
        fa = fasta.Fasta(path, 'w')
        fa.write(fasta.FastaEntry('>header', 'ACTG'))
        fa.close()
        with open(path, 'r') as f:
            self.assertEqual('>header\nACTG\n', f.read())
