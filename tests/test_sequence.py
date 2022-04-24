from collections import Counter
from unittest import TestCase, mock

import numpy as np

from ngs_tools import sequence

from . import mixins


class TestSequence(mixins.TestMixin, TestCase):

    def test_alignment_to_cigar(self):
        self.assertEqual('4D', sequence.alignment_to_cigar('ACGT', '----'))
        self.assertEqual('1M2D1M', sequence.alignment_to_cigar('ACGT', 'A--T'))
        self.assertEqual('4I', sequence.alignment_to_cigar('----', 'AAAT'))
        self.assertEqual(
            '3M1X', sequence.alignment_to_cigar('ACGT', 'ACNV', mismatch=True)
        )
        with self.assertRaises(sequence.SequenceError):
            sequence.alignment_to_cigar('AAAA', 'AAA')
        with self.assertRaises(sequence.SequenceError):
            sequence.alignment_to_cigar('----', '----')

    def test_call_consensus_with_qualities(self):
        with open(self.sequences_path, 'r') as f1, open(self.qualities_path,
                                                        'r') as f2:
            sequences = [line.strip() for line in f1 if not line.isspace()]
            qualities = [line.strip() for line in f2 if not line.isspace()]

        consensuses, assignments = sequence.call_consensus_with_qualities(
            sequences, qualities
        )
        counts = Counter(assignments)

        # Check that assignments are ordered
        self.assertEqual(
            sorted(counts.values(), reverse=True),
            [counts[i] for i in range(len(counts))]
        )

        common = counts.most_common(2)
        self.assertTrue(common[0][1] > 50)
        self.assertTrue(common[1][1] > 25)

    def test_call_consensus_with_qualities_allow_ambiguous(self):
        sequences = ['AAAC', 'AAAG']
        qualities = ['AAAA', 'AAAA']

        consensuses, assignments = sequence.call_consensus_with_qualities(
            sequences, qualities, allow_ambiguous=True
        )
        self.assertEqual(consensuses, ['AAAS'])
        np.testing.assert_equal(assignments, [0, 0])

    def test_call_consensus_with_qualities_return_qualities(self):
        with open(self.sequences_path, 'r') as f1, open(self.qualities_path,
                                                        'r') as f2:
            sequences = [line.strip() for line in f1 if not line.isspace()]
            qualities = [line.strip() for line in f2 if not line.isspace()]

        consensuses, assignments, consensus_qualities = sequence.call_consensus_with_qualities(
            sequences, qualities, return_qualities=True
        )
        counts = Counter(assignments)
        common = counts.most_common(2)
        self.assertTrue(common[0][1] > 50)
        self.assertTrue(common[1][1] > 25)

    def test_levenshtein_distance(self):
        self.assertEqual(1, sequence.levenshtein_distance('AC', 'AT'))
        self.assertEqual(0, sequence.levenshtein_distance('AT', 'AN'))
        self.assertEqual(2, sequence.levenshtein_distance('XZ', 'ZX'))

    def test_levenshtein_distance_raises_error(self):
        with mock.patch('ngs_tools.sequence.LEVENSHTEIN_DISTANCE_ALIGNER',
                        None):
            with self.assertRaises(sequence.SequenceError):
                sequence.levenshtein_distance('AC', 'AT')

    def test_hamming_distance(self):
        self.assertEqual(0, sequence.hamming_distance('ACTG', 'ACTG'))
        self.assertEqual(1, sequence.hamming_distance('ACTG', 'ACTT'))
        self.assertEqual(0, sequence.hamming_distance('ACTG', 'ACTN'))

    def test_hamming_distances(self):
        np.testing.assert_equal(
            np.array([0, 1, 0]),
            sequence.hamming_distances('ACTG', ['ACTG', 'ACTT', 'ACTN'])
        )

    def test_hamming_distance_matrix(self):
        np.testing.assert_equal(
            np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
            sequence.hamming_distance_matrix(['ACTG', 'ACTT', 'ACTN'],
                                             ['ACTG', 'ACTT', 'ACTN'])
        )

    def test_pairwise_hamming_distances(self):
        np.testing.assert_equal(
            np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
            sequence.pairwise_hamming_distances(['ACTG', 'ACTT', 'ACTN'])
        )

    def test_correct_sequences_to_whitelist(self):
        sequences = ['ACTG', 'ACTT', 'AGCC', 'TTTT']
        qualities = ['AAAA', 'AAAA', 'AAAA', 'AAAA']
        whitelist = ['ACTG', 'TTTN']
        with mock.patch('ngs_tools.sequence.utils.progress', mixins.tqdm_mock),\
            mock.patch('ngs_tools.sequence.progress', mixins.tqdm_mock):
            corrections = sequence.correct_sequences_to_whitelist(
                sequences, qualities, whitelist
            )
        self.assertEqual(['ACTG', 'ACTG', None, 'TTTN'], corrections)

    def test_correct_sequences_to_whitelist_simple(self):
        sequences = ['ACTG', 'ACTT', 'AGCC', 'TTTT']
        whitelist = ['ACTG', 'TTTN']
        with mock.patch('ngs_tools.sequence.utils.progress', mixins.tqdm_mock),\
            mock.patch('ngs_tools.sequence.progress', mixins.tqdm_mock):
            corrections = sequence.correct_sequences_to_whitelist_simple(
                sequences, whitelist
            )
        self.assertEqual({
            'ACTG': 'ACTG',
            'ACTT': 'ACTG',
            'AGCC': None,
            'TTTT': 'TTTN'
        }, corrections)
