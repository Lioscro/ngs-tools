import array
from typing import List, Union

import numpy as np
import pysam

NUCLEOTIDES_STRICT = ['A', 'C', 'G', 'T']
NUCLEOTIDES_PERMISSIVE = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
NUCLEOTIDES = NUCLEOTIDES_STRICT + NUCLEOTIDES_PERMISSIVE
NUCLEOTIDES_AMBIGUOUS = {
    'N': ('A', 'C', 'G', 'T'),
    'R': ('A', 'G'),
    'Y': ('C', 'T'),
    'S': ('G', 'C'),
    'W': ('A', 'T'),
    'K': ('G', 'T'),
    'M': ('A', 'C'),
    'B': ('C', 'G', 'T'),
    'D': ('A', 'G', 'T'),
    'H': ('A', 'C', 'T'),
    'V': ('A', 'C', 'G'),
}
NUCLEOTIDE_COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N',
    'R': 'Y',
    'Y': 'R',
    'S': 'W',
    'W': 'S',
    'K': 'M',
    'M': 'K',
    'B': 'V',
    'D': 'H',
    'H': 'D',
    'V': 'B',
}


def _sequence_to_array(sequence, chars, l=None):  # noqa: E741
    char_to_idx = {c: i for i, c in enumerate(chars)}
    arr = np.zeros((len(chars), l or len(sequence)), dtype=bool)
    for i, c in enumerate(sequence):
        arr[char_to_idx[c], i] = True
    return arr


def _most_likely_sequence(positional_probs, chars):
    indices = positional_probs.argmax(axis=0)
    return ''.join(chars[i] for i in indices)


def _calculate_positional_probs(sequences, qualities):
    positional_probs = np.zeros(sequences[0].shape)
    for seq, qual in zip(sequences, qualities):
        np.add(positional_probs, qual, out=positional_probs, where=seq)
    return positional_probs


def call_consensus_with_qualities(
    sequences: List[str],
    qualities: Union[List[str], List[array.array]],
    q_threshold: int = 30,
    proportion: float = 0.05,
    return_qualities: bool = False,
):
    """Given a list of sequences and their base qualities, constructs a *set* of consensus
    sequences by iteratively constructing a consensus (by selecting the most likely
    base at each position) and assigning sequences with match probability <=
    max(min(match probability), `q_threshold` * (`proportion` * length of longest sequence))
    to this consensus. Then, the consensus is updated by constructing the consensus only
    among these sequences. The match probability of a sequence to a consensus is the sum of
    the quality values where they do not match (equivalent to negative log probability that
    all mismatches were sequencing errors). Provided work well for most cases.
    """
    # Check number of sequences and their lengths match with provided qualities
    if len(sequences) != len(qualities):
        raise Exception(
            f'{len(sequences)} sequences and {len(qualities)} qualities were provided'
        )
    if any(len(seq) != len(qual) for seq, qual in zip(sequences, qualities)):
        raise Exception(
            'length of each sequence must match length of each quality string'
        )

    def _call_consensus(seqs, quals, chars, thresh):
        if len(seqs) == 1:
            return _most_likely_sequence(seqs[0], chars), np.array([True],
                                                                   dtype=bool)

        positional_probs = _calculate_positional_probs(seqs, quals)
        consensus_indices = positional_probs.argmax(axis=0)

        # For each sequence, calculate the probability that the sequence was actually
        # equal the consensus, but the different bases are due to sequencing errors
        # NOTE: should we also be looking at probability that matches are correct?
        probs = []
        for seq, qual in zip(seqs, quals):
            p = np.sum(qual[consensus_indices != seq.argmax(axis=0)])
            probs.append(p)
        probs = np.array(probs)
        assigned = probs <= max(thresh, min(probs))

        # NOTE: we construct a new consensus from assigned sequences
        assigned_seqs = seqs[assigned]
        assigned_quals = quals[assigned]
        return _most_likely_sequence(
            _calculate_positional_probs(assigned_seqs, assigned_quals), chars
        ), assigned

    # Convert sequences to array representations
    chars = sorted(set.union(*[set(s) for s in sequences]))
    l = max(len(s) for s in sequences)  # noqa: E741
    sequences_arrays = np.array([
        _sequence_to_array(sequence, chars, l=l) for sequence in sequences
    ])
    # Convert quality strings to quality values (integers)
    qualities_arrays = []
    for quals in qualities:
        arr = np.array(
            pysam.
            qualitystring_to_array(quals) if isinstance(quals, str) else quals
        )
        arr.resize(l)
        qualities_arrays.append(arr)
    qualities_arrays = np.array(qualities_arrays)

    # Iteratively call consensus sequences. This used to be done recursively, but there were cases
    # when Python's recursion limit would be reached. Thankfully, all recursive algorithms can be
    # rewritten to be iterative.
    threshold = q_threshold * (l * proportion)
    consensuses = []
    assignments = np.full(len(sequences), -1, dtype=int)
    index_transform = {i: i for i in range(len(sequences))}
    _sequences_arrays = sequences_arrays.copy()
    _qualities_arrays = qualities_arrays.copy()
    while True:
        consensus, assigned = _call_consensus(
            _sequences_arrays, _qualities_arrays, chars, threshold
        )
        if consensus in consensuses:
            label = consensuses.index(consensus)
        else:
            label = len(consensuses)
            consensuses.append(consensus)

        assigned_indices = assigned.nonzero()[0]
        unassigned_indices = (~assigned).nonzero()[0]
        assignments[[index_transform[i] for i in assigned_indices]] = label
        if all(assigned):
            break
        index_transform = {
            i: index_transform[j]
            for i, j in enumerate(unassigned_indices)
        }
        _sequences_arrays = _sequences_arrays[~assigned]
        _qualities_arrays = _qualities_arrays[~assigned]

    # Compute qualities for each consensus sequence if return_qualities = True
    if return_qualities:
        consensuses_qualities = []
        for i, consensus in enumerate(consensuses):
            assigned = assignments == i
            assigned_sequences = sequences_arrays[assigned]
            assigned_qualities = qualities_arrays[assigned]
            consensus_array = _sequence_to_array(consensus, chars, l)

            # (assigned_sequences & consensus_array) is a 3-dimensional array
            # First dimension contains each sequence, second contains base identity,
            # third contains positions, so (assigned_sequences & consensus_array) contains True
            # in positions of each sequence where the sequence has the same base as the
            # consensus. Taking the any(axis=1) of this gives a 2D matrix where each row
            # corresponds to each sequence and each column contains True if the base at that
            # position in that sequence matches the consensus. Multiplying this
            # boolean mask with the assigned_qualities gives the base quality of each
            # sequence only at the positions where the base matches that of the consensus.
            # Then, we take the maximum quality among all these bases.
            consensus_qualities = (
                assigned_qualities *
                (assigned_sequences & consensus_array).any(axis=1)
            ).max(axis=0)
            consensuses_qualities.append(consensus_qualities)
        return consensuses, assignments, consensus_qualities
    else:
        return consensuses, assignments
