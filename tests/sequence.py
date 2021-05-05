from collections import Counter
from unittest import mock, TestCase

import numpy as np

from ngsutils import sequence

from tests.mixins import TestMixin, dummy_function


class TestSequence(TestMixin, TestCase):

    def test_call_consensus_with_qualities(self):
        with open(self.sequences_path, 'r') as f1, open(self.qualities_path,
                                                        'r') as f2:
            sequences = [line.strip() for line in f1 if not line.isspace()]
            qualities = [line.strip() for line in f2 if not line.isspace()]

        consensuses, assignments = sequence.call_consensus_with_qualities(
            sequences, qualities
        )
        counts = Counter(assignments)
        common = counts.most_common(2)
        self.assertTrue(common[0][1] > 50)
        self.assertTrue(common[1][1] > 25)

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
        with mock.patch('ngsutils.sequence.utils.tqdm', dummy_function),\
            mock.patch('ngsutils.sequence.tqdm', dummy_function):
            corrections = sequence.correct_sequences_to_whitelist(
                sequences, qualities, whitelist
            )
        self.assertEqual(['ACTG', 'ACTG', None, 'TTTN'], corrections)
