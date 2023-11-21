import array
import re
from collections import Counter
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pysam
from joblib import delayed
from numba import njit

try:
    from pyseq_align import NeedlemanWunsch
except ImportError:
    NeedlemanWunsch = None

from . import utils
from .logging import logger
from .progress import progress

NUCLEOTIDES_STRICT = ['A', 'C', 'G', 'T']
NUCLEOTIDES_PERMISSIVE = [
    'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '-'
]
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
    '-': tuple(),
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
NUCLEOTIDE_MASKS = {
    n: np.array([
        _n in NUCLEOTIDES_AMBIGUOUS.get(n, [n]) for _n in NUCLEOTIDES_STRICT
    ],
                dtype=bool)
    for n in NUCLEOTIDES
}
MASK_TO_NUCLEOTIDE = {tuple(mask): n for n, mask in NUCLEOTIDE_MASKS.items()}
LEVENSHTEIN_DISTANCE_ALIGNER = NeedlemanWunsch(
    match=0,
    mismatch=-1,
    gap_open=0,
    gap_extend=-1,
    substitution_matrix={
        n: {
            _n: int((NUCLEOTIDE_MASKS[n] & NUCLEOTIDE_MASKS[_n]).any()) - 1
            for _n in NUCLEOTIDES
        }
        for n in NUCLEOTIDES
    }
) if NeedlemanWunsch is not None else None
SEQUENCE_PARSER = re.compile(r'[^atcgATCG]')


class SequenceError(Exception):
    pass


def alignment_to_cigar(
    reference: str, query: str, mismatch: bool = False
) -> str:
    """Convert an alignment to a CIGAR string.

    The CIGAR is always constructed relative to the `reference` (i.e. as
    insertions/deletions from the reference). The provided sequences must
    represent their alignments, containing the "-" character at positions where
    there are gaps in one or another.

    Args:
        reference: The reference sequence alignment
        query: The query sequence alignment
        mismatch: Whether or not to use the "X" CIGAR operation in places of
            mismatches. Defaults to `True`, such that all non-gaps are considered
            matches with the "M" CIGAR operation. When there are ambiguous
            characters, a mismatch occurs only when the two sets of possible
            nucleotides are exclusive.

    Returns:
        The CIGAR string representing the provided alignment

    Raises:
        SequenceError: if the `reference` and `query` do not have the same
            length, or if there are gaps in both alignments at the same position
    """
    if len(reference) != len(query):
        raise SequenceError('The two alignments must have the same lengths.')

    cigar = ''
    prev_state = None
    count = 0
    for r, q in zip(reference, query):
        if r == '-' and q == '-':
            raise SequenceError(
                'The two alignments have gaps at the same position.'
            )

        if r == '-':
            state = 'I'
        elif q == '-':
            state = 'D'
        elif mismatch:
            _r = set(NUCLEOTIDES_AMBIGUOUS.get(r, (r,)))
            _q = set(NUCLEOTIDES_AMBIGUOUS.get(q, (q,)))
            if _r & _q:
                state = 'M'
            else:
                state = 'X'
        else:
            state = 'M'

        if state == prev_state:
            count += 1
        else:
            if prev_state is not None:
                cigar += f'{count}{prev_state}'
            prev_state = state
            count = 1

    if prev_state is not None:
        cigar += f'{count}{prev_state}'
    return cigar


def complement_sequence(sequence: str, reverse: bool = False) -> str:
    """Complement the given sequence, with optional reversing.

    Args:
        sequence: Input sequence
        reverse: Whether or not to perform reverse complementation

    Returns:
        Complemented (and optionally reversed) string
    """
    sequence = sequence.upper()
    if reverse:
        sequence = reversed(sequence)
    return ''.join(NUCLEOTIDE_COMPLEMENT[c] for c in sequence)


def _sequence_to_array(
    sequence: str,
    l: Optional[int] = None  # noqa: E741
) -> np.ndarray:  # noqa: E741
    sequence = sequence.upper()
    for c in sequence:
        if c not in NUCLEOTIDES:
            raise SequenceError(f'Unknown nucleotide `{c}`')

    arr = np.zeros((len(NUCLEOTIDES_STRICT), l or len(sequence)), dtype=bool)
    for i, c in enumerate(sequence):
        arr[NUCLEOTIDE_MASKS[c], i] = True
    return arr


def _qualities_to_array(
    qualities: Union[str, array.array],
    l: Optional[int] = None  # noqa: E741
) -> np.ndarray:
    if l and l < len(qualities):
        raise SequenceError('`l` can not be smaller than length of `qualities`')

    arr = np.array(
        pysam.qualitystring_to_array(qualities)
        if isinstance(qualities, str) else qualities,
        dtype=np.uint8
    )
    if l:
        arr.resize(l)
    return arr


@njit
def _most_likely_array(positional_probs: np.ndarray) -> np.ndarray:
    most_likely = np.zeros(positional_probs.shape, dtype=np.bool_)
    for i in range(positional_probs.shape[1]):
        most_likely[:, i] = positional_probs[:, i] == positional_probs[:,
                                                                       i].max()
    return most_likely


def _most_likely_sequence(
    positional_probs: np.ndarray, allow_ambiguous: bool = False
) -> str:
    if not allow_ambiguous:
        indices = positional_probs.argmax(axis=0)
        return ''.join(NUCLEOTIDES_STRICT[i] for i in indices)
    else:
        arr = _most_likely_array(positional_probs)
        s = ''
        for i in range(arr.shape[1]):
            s += MASK_TO_NUCLEOTIDE[tuple(arr[:, i])]
        return s


def _disambiguate_sequence(sequence: np.ndarray) -> List[str]:
    sequences = ['']
    for pos in sequence.T:
        new_sequences = []
        for i in pos.nonzero()[0]:
            new_sequences += [
                f'{seq}{NUCLEOTIDES_STRICT[i]}' for seq in sequences
            ]
        sequences = new_sequences
    return sequences


@njit
def _calculate_positional_probs(
    sequences: np.ndarray, qualities: np.ndarray
) -> np.ndarray:
    positional_probs = np.zeros(sequences[0].shape, dtype=np.uint)
    for seq, qual in zip(sequences, qualities):
        _qual = qual.repeat(seq.shape[0]).reshape(-1, seq.shape[0]).T
        positional_probs += _qual * seq
    return positional_probs


def call_consensus(sequences: List[str], proportion: float = 0.05):
    """Call consensus sequences from a set of sequences. Internally, this
    function calls :func:`call_consensus_with_qualities` with all qualities set
    to 1 and ``q_threshold`` set to 1. See the documentation of this function for
    details on how consensuses are called.

    Args:
        sequences: Sequences to call consensus for
        proportion: Proportion of each sequence to allow mismatched bases to be
            above ``q_threshold``

    Returns:
        List of consensus sequences
        Numpy array of assignments for each sequence in ``sequences``
    """
    qualities = [np.full(len(seq), 1) for seq in sequences]
    return call_consensus_with_qualities(
        sequences,
        qualities,
        q_threshold=1,
        proportion=proportion,
        return_qualities=False
    )


def call_consensus_with_qualities(
    sequences: List[str],
    qualities: Union[List[str], List[array.array], List[np.ndarray]],
    q_threshold: Optional[int] = None,
    proportion: float = 0.05,
    allow_ambiguous: bool = False,
    return_qualities: bool = False,
) -> Union[Tuple[List[str], np.ndarray], Tuple[List[str], np.ndarray,
                                               List[str]]]:
    """Given a list of sequences and their base qualities, constructs a *set* of consensus
    sequences by iteratively constructing a consensus (by selecting the most likely
    base at each position) and assigning sequences with match probability <=
    max(min(match probability), `q_threshold` * (`proportion` * length of longest sequence))
    to this consensus. Then, the consensus is updated by constructing the consensus only
    among these sequences. The match probability of a sequence to a consensus is the sum of
    the quality values where they do not match (equivalent to negative log probability that
    all mismatches were sequencing errors).

    Note:
        This function does not perform any alignment among consensus sequences.
        To detect any insertions/deletions, call this function and then
        perform alignment among the called consensus sequences.

    Args:
        sequences: Sequences to call consensus for
        qualities: Quality scores for the sequences
        q_threshold: Quality threshold. Defaults to the median quality score among
            all bases in all sequences.
        proportion: Proportion of each sequence to allow mismatched bases to be
            above ``q_threshold``. Defaults to 0.05.
        allow_ambiguous: Allow ambiguous bases in the consensus sequences. Defaults to
            False, which, on ties, selects a single base lexicographically. This
            option only has an effect when constructing the final consensus
            sequence as a string, not when calculating error probabilities.
        return_qualities: Whether or not to return qualities for the consensuses.
            Defaults to False.

    Returns:
        List of consensus sequences
        Numpy array of assignments for each sequence in ``sequences``
        Qualities for each of the consensus sequences, if ``return_qualities`` is True

    Raises:
        SequenceError: if any sequence-quality pair have different lengths or
            number of provided sequences does not match number of provided
            qualities
    """
    # Check number of sequences and their lengths match with provided qualities
    if len(sequences) != len(qualities):
        raise SequenceError(
            f'{len(sequences)} sequences and {len(qualities)} qualities were provided'
        )
    if any(len(seq) != len(qual) for seq, qual in zip(sequences, qualities)):
        raise SequenceError(
            'length of each sequence must match length of each quality string'
        )

    def _consensus_probs(consensus_array, seqs, quals):
        probs = np.zeros(seqs.shape[0])
        for i in range(seqs.shape[0]):
            probs[i] = quals[i][_mismatch_mask(consensus_array, seqs[i])].sum()
        return probs

    _consensus_probs_njit = njit(_consensus_probs)

    def _call_consensus(seqs, quals, thresh):
        if len(seqs) == 1:
            return seqs[0], np.array([0],
                                     dtype=np.uint), np.array([], dtype=np.uint)

        positional_probs = _calculate_positional_probs(seqs, quals)
        consensus_array = _most_likely_array(positional_probs)

        # Numbaized function is efficient when there are many sequences
        probs_func = _consensus_probs if len(
            seqs
        ) < 100 else _consensus_probs_njit

        # For each sequence, calculate the probability that the sequence was actually
        # equal the consensus, but the different bases are due to sequencing errors
        # NOTE: should we also be looking at probability that matches are correct?
        probs = probs_func(consensus_array, seqs, quals)
        # the max is taken to deal with case where there is only one sequence
        assigned = probs <= max(thresh, min(probs))
        assigned_indices = assigned.nonzero()[0]
        unassigned_indices = (~assigned).nonzero()[0]

        # NOTE: we construct a new consensus from assigned sequences
        assigned_seqs = seqs[assigned_indices]
        assigned_quals = quals[assigned_indices]
        return _calculate_positional_probs(
            assigned_seqs, assigned_quals
        ), assigned_indices, unassigned_indices

    if len(sequences) > 10000:
        logger.warning(
            'More than 10000 sequences provided. This may take a while.'
        )

    # Convert sequences to array representations
    l = max(len(s) for s in sequences)  # noqa: E741
    sequences_arrays = np.array([
        _sequence_to_array(sequence, l=l) for sequence in sequences
    ])
    # Convert quality strings to quality values (integers)
    qualities_arrays = np.array([
        _qualities_to_array(quals, l=l) for quals in qualities
    ])
    # Compute dynamic quality threshold if not provided
    if not q_threshold:
        q_threshold = np.median(qualities_arrays)

    # Iteratively call consensus sequences. This used to be done recursively, but there were cases
    # when Python's recursion limit would be reached. Thankfully, all recursive algorithms can be
    # rewritten to be iterative.
    threshold = q_threshold * (l * proportion)
    consensus_index = {}
    assignments = np.full(len(sequences), -1, dtype=int)
    index_transform = {i: i for i in range(len(sequences))}
    _sequences_arrays = sequences_arrays.copy()
    _qualities_arrays = qualities_arrays.copy()
    while True:
        probs, assigned_indices, unassigned_indices = _call_consensus(
            _sequences_arrays, _qualities_arrays, threshold
        )

        consensus = _most_likely_sequence(probs, allow_ambiguous)
        label = consensus_index.setdefault(consensus, len(consensus_index))

        assignments[[index_transform[i] for i in assigned_indices]] = label
        if unassigned_indices.size == 0:
            break
        index_transform = {
            i: index_transform[j]
            for i, j in enumerate(unassigned_indices)
        }
        _sequences_arrays = _sequences_arrays[unassigned_indices]
        _qualities_arrays = _qualities_arrays[unassigned_indices]

    # Reorder assignment so that consensuses[0] is the consensus with the greatest
    # number of assigned sequences, and so on.
    consensuses = list(consensus_index.keys())
    reordered_consensuses = []
    reordered_assignments = np.full(len(sequences), -1, dtype=int)
    for i, _ in Counter(assignments).most_common():
        reordered_consensuses.append(consensuses[i])
        reordered_assignments[assignments == i] = len(reordered_consensuses) - 1
    consensuses = reordered_consensuses
    assignments = reordered_assignments

    # Compute qualities for each consensus sequence if return_qualities = True
    if return_qualities:
        consensuses_qualities = []
        for i, consensus in enumerate(consensuses):
            assigned = assignments == i
            assigned_sequences = sequences_arrays[assigned]
            assigned_qualities = qualities_arrays[assigned]
            consensus_array = _sequence_to_array(consensus, l)

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
            consensuses_qualities.append(
                pysam.qualities_to_qualitystring(consensus_qualities)
            )
        return consensuses, assignments, consensuses_qualities
    else:
        return consensuses, assignments


def levenshtein_distance(sequence1: str, sequence2: str) -> int:
    """Calculate the Levenshtein (edit) distance between two sequences.
    This is calculated by calling :func:`pyseq_align.NeedlemanWunsch.align`
    with the appropriate scores/penalties.

    Args:
        sequence1: First sequence
        sequence2: Second sequence

    Returns:
        Levenshtein distance
    """
    if LEVENSHTEIN_DISTANCE_ALIGNER is None:
        raise SequenceError(
            'Levenshtein distance requires the `pyseq-align` library, which '
            'can be installed with `pip install pyseq-align` or '
            '`pip install ngs-tools[levenshtein]`.'
        )
    return -LEVENSHTEIN_DISTANCE_ALIGNER.align(sequence1, sequence2).score


@njit
def _mismatch_mask(sequence1: np.ndarray, sequence2: np.ndarray) -> int:
    not_and = ~(sequence1 & sequence2)
    result = np.ones(not_and.shape[1], dtype=np.bool_)
    for i in range(not_and.shape[0]
                   ):  # np.all with arguments isn't available in @njit
        result &= not_and[i]
    return result


@njit
def _mismatch_masks(
    sequence: np.ndarray, whitelist: np.ndarray, d: int
) -> np.ndarray:
    indices = []
    masks = []
    for i, bc_array in enumerate(whitelist):
        mask = _mismatch_mask(sequence, bc_array)
        if mask.sum() > d:
            continue
        indices.append(i)
        masks.append(mask)
    return indices, masks


@njit
def _hamming_distance(sequence1: np.ndarray, sequence2: np.ndarray) -> int:
    return _mismatch_mask(sequence1, sequence2).sum()


def hamming_distance(sequence1: str, sequence2: str) -> int:
    """Calculate the hamming distance between two sequences.

    Args:
        sequence1: First sequence
        sequence2: Second sequence

    Returns:
        Hamming distance

    Raises:
        SequenceError: When the sequences are unequal lengths
    """
    if len(sequence1) != len(sequence2):
        raise SequenceError('Unequal lengths')

    return _hamming_distance(
        _sequence_to_array(sequence1), _sequence_to_array(sequence2)
    )


@njit
def _hamming_distances(
    sequence: np.ndarray, sequences: np.ndarray
) -> np.ndarray:
    distances = np.zeros(sequences.shape[0], dtype=np.uint)
    for i, seq in enumerate(sequences):
        distances[i] = _hamming_distance(sequence, seq)
    return distances


def hamming_distances(sequence: str, sequences: List[str]) -> np.ndarray:
    """Calculate the hamming distance between a sequence and a list of sequences.

    Args:
        sequence: Sequence
        sequences: List of sequences

    Returns:
        Numpy array of hamming distances, where each index ``i`` corresponds to
        the hamming distance between ``sequence`` and ``sequences[i]``

    Raises:
        SequenceError: When any of the sequences are unequal length
    """
    if any(len(sequence) != len(seq) for seq in sequences):
        raise SequenceError('All sequences must be equal length')

    sequence = _sequence_to_array(sequence)
    sequences = np.array([_sequence_to_array(seq) for seq in sequences])
    return _hamming_distances(sequence, sequences)


@njit
def _hamming_distance_matrix(
    sequences1: np.ndarray, sequences2: np.ndarray
) -> np.ndarray:
    distances = np.zeros((sequences1.shape[0], sequences2.shape[0]),
                         dtype=np.uint)
    for i, seq1 in enumerate(sequences1):
        distances[i] = _hamming_distances(seq1, sequences2)
    return distances


def hamming_distance_matrix(
    sequences1: List[str], sequences2: List[str]
) -> np.ndarray:
    """Calculate all pairwise hamming distances between two lists of sequences.

    Args:
        sequences1: List of sequences
        sequences2: List of sequences

    Returns:
        Numpy array of hamming distances, where each index ``(i, j)`` corresponds to
        the hamming distance between ``sequences1[i]`` and ``sequences2[j]``

    Raises:
        SequenceError: When any of the sequences are unequal length
    """
    if any(len(sequences1[0]) != len(seq) for seq in sequences1 + sequences2):
        raise SequenceError('All sequences must be equal length')

    sequences1 = np.array([
        _sequence_to_array(sequence) for sequence in sequences1
    ])
    sequences2 = np.array([
        _sequence_to_array(sequence) for sequence in sequences2
    ])
    return _hamming_distance_matrix(sequences1, sequences2)


@njit
def _pairwise_hamming_distances(sequences: np.ndarray) -> np.ndarray:
    distances = np.zeros((len(sequences), len(sequences)), dtype=np.uint)
    for i in range(sequences.shape[0]):
        for j in range(i, sequences.shape[0]):
            d = _hamming_distance(sequences[i], sequences[j])
            distances[i, j] = d
            distances[j, i] = d
    return distances


def pairwise_hamming_distances(sequences: List[str]) -> np.ndarray:
    """Calculate all pairwise hamming distances between combinations of sequences
    from a single list.

    Args:
        sequences: List of sequences

    Returns:
        Numpy array of hamming distances, where each index ``(i, j)`` corresponds to
        the hamming distance between ``sequences[i]`` and ``sequences[j]``

    Raises:
        SequenceError: When any of the sequences are unequal length
    """
    if any(len(sequences[0]) != len(seq) for seq in sequences):
        raise SequenceError('All sequences must be equal length')
    sequences = np.array([
        _sequence_to_array(sequence) for sequence in sequences
    ])
    return _pairwise_hamming_distances(sequences)


@njit
def _correct_to_whitelist(
    qualities: np.ndarray,
    indices: np.ndarray,
    masks: np.ndarray,
    log10_proportions: np.ndarray,
) -> Tuple[int, float]:
    best_bc = -1
    max_log10_likelihood = -np.inf
    log10_likelihoods = []
    for i, mask in zip(indices, masks):
        log10p_edit = -(qualities[mask].sum() / 10)
        log10_likelihood = log10_proportions[i] + log10p_edit
        log10_likelihoods.append(log10_likelihood)

        if log10_likelihood > max_log10_likelihood:
            max_log10_likelihood = log10_likelihood
            best_bc = i

    log10_confidence = 0
    if best_bc >= 0:
        log10_confidence = max_log10_likelihood - np.log10(
            (np.power(10, np.array(log10_likelihoods))).sum()
        )
    return best_bc, log10_confidence


def correct_sequences_to_whitelist(
    sequences: List[str],
    qualities: Union[List[str], List[array.array]],
    whitelist: List[str],
    d: int = 1,
    confidence: float = 0.95,
    n_threads: int = 1,
    show_progress: bool = False,
) -> List[Union[str, None]]:
    """Correct a list of sequences to a whitelist within `d` hamming distance.
    Note that `sequences` can contain duplicates, but `whitelist` can not.

    For a given sequence, if there are multiple barcodes in the whitelist to which
    its distance is <= `d`, the sequence is assigned to a barcode by using the
    prior probability that the sequence originated from the barcode. If the confidence
    that the sequence originated from the most likely barcode is less than `confidence`,
    assignment is skipped.

    This procedure follows the barcode correction procedure in Cell Ranger by
    10X Genomics. Some modifications were made to support ambiguous bases and
    prevent floating-point underflow.
    https://kb.10xgenomics.com/hc/en-us/articles/115003822406-How-does-Cell-Ranger-correct-barcode-sequencing-errors

    Note:
        Only hamming distance is supported (not Levenshtein distance).

    Args:
        sequences: List of sequence strings
        qualities: List of quality strings or list of array of qualities (as
            returned by :func:`pysam.qualitystring_to_array`)
        whitelist: List of whitelist sequences to correct to
        d: Hamming distance threshold. Sequences will be corrected to the whitelist
            with hamming distance <= ``d``. Defaults to 1.
        confidence: When a sequence has the same minimum hamming distance to
            multiple whitelisted sequences, the sequence is assigned to the
            best whitelisted sequence (using prior probabilities) if the likelihood
            ratio of this and the sum of all likelihoods is greater than or equal to
            this value. Defaults to 0.95.
        n_threads: Number of threads to use. Defaults to 1.
        show_progress: Whether to display a progress bar. Defaults to True.

    Raises:
        SequenceError: If all the lengths of each sequence, qualitiy and
            whitelisted sequence are not equal, the number of sequences and
            qualities provided are not equal or the whitelist contains duplicates.

    Returns:
        The corrections as a list of whitelisted sequences. For sequences that
        could not be corrected, the corresponding position contains None.
    """
    # Check number of sequences and their lengths match with provided qualities
    if len(sequences) != len(qualities):
        raise SequenceError(
            f'{len(sequences)} sequences and {len(qualities)} qualities were provided'
        )
    if any(len(seq) != len(qual) for seq, qual in zip(sequences, qualities)):
        raise SequenceError(
            'length of each sequence must match length of each quality string'
        )
    if len(set(whitelist)) != len(whitelist):
        raise SequenceError('`whitelist` contains duplicates')
    for seq in sequences:
        if len(seq) != len(sequences[0]):
            raise SequenceError(
                'all sequences in `sequences` must be of same length'
            )
    for bc in whitelist:
        if len(bc) != len(whitelist[0]):
            raise SequenceError(
                'all sequences in `whitelist` must be of same length'
            )
    if len(sequences[0]) != len(whitelist[0]):
        raise SequenceError(
            'all sequences in `sequences` and `whitelist` must be of same length'
        )

    # Step 1: any exact matches (ignoring ambiguities) can be assigned Immediately
    counts = Counter(sequences)
    whitelist_counts = np.zeros(len(whitelist), dtype=float)
    whitelist_indices = {bc: i for i, bc in enumerate(whitelist)}
    unmatched_sequences = []
    matches = {}
    for seq in progress(list(counts.keys()), desc='[1/4] Finding exact matches',
                        disable=not show_progress):
        if seq in whitelist_indices:
            matches[seq] = seq
            whitelist_counts[whitelist_indices[seq]] += counts[seq]
        else:
            unmatched_sequences.append(seq)

    # Step 2: Construct whitelist mask
    whitelist_arrays = np.zeros(
        (len(whitelist), len(NUCLEOTIDES_STRICT), len(whitelist[0])),
        dtype=bool
    )
    for i, bc in enumerate(progress(whitelist, desc='[2/4] Constructing masks',
                                    disable=not show_progress)):
        whitelist_arrays[i] = _sequence_to_array(bc)

    # Step 3: Find all mismatch masks for which hamming distance <= d
    mismatch_cache = {}
    for i, (indices, masks) in enumerate(utils.ParallelWithProgress(
            n_jobs=n_threads, total=len(unmatched_sequences),
            desc='[3/4] Finding mismatches', disable=not show_progress
    )(delayed(_mismatch_masks)(_sequence_to_array(seq), whitelist_arrays, d=d)
      for seq in unmatched_sequences)):
        indices = np.array(indices, dtype=int)
        masks = np.array(masks, dtype=bool)
        seq = unmatched_sequences[i]
        mismatch_cache[seq] = (indices, masks)

        if indices.shape[0] == 0:
            continue

        # Check if there is an exact match to more than one barcode. If there was,
        # don't add the count to the whitelisted count.
        match_indices = indices[masks.sum(axis=1) == 0]
        if len(match_indices) == 1:
            matches[seq] = whitelist[match_indices[0]]
        if len(match_indices) > 0:
            whitelist_counts[match_indices] += counts[seq] / len(match_indices)

    # Step 4: correct all other sequences to whitelist
    corrections = [None] * len(sequences)
    with progress(total=len(sequences), desc='[4/4] Correcting sequences',
                  disable=not show_progress) as pbar:
        for i, sequence in enumerate(sequences):
            if sequence in matches:
                corrections[i] = matches[sequence]
                pbar.update(1)

        # Calculate proportions
        whitelist_pseudo = sum(whitelist_counts) + len(whitelist)
        whitelist_log10_proportions = np.log10(
            (whitelist_counts + 1) / whitelist_pseudo
        )

        confidence = np.log10(confidence)
        for i, seq in enumerate(sequences):
            if corrections[i] is not None:
                continue

            best_bc, log10_confidence = _correct_to_whitelist(
                _qualities_to_array(qualities[i]), mismatch_cache[seq][0],
                mismatch_cache[seq][1].reshape(1,
                                               -1), whitelist_log10_proportions
            )
            if best_bc >= 0 and log10_confidence >= confidence:
                corrections[i] = whitelist[best_bc]
            pbar.update(1)
            pbar.refresh()

    return corrections


def correct_sequences_to_whitelist_simple(
    sequences: List[str],
    whitelist: List[str],
    d: int = 1,
    n_threads: int = 1,
    show_progress: bool = False,
) -> Dict[str, Union[str, None]]:
    """Correct a list of sequences to a whitelist within `d` hamming distance.
    Note that `sequences` can contain duplicates, but `whitelist` can not.
    Unlike :func:`correct_sequences_to_whitelist`, this function takes a naive
    approach and discards any sequences that can be corrected to multiple
    whitelisted sequences.

    Args:
        sequences: List of sequence strings
        whitelist: List of whitelist sequences to correct to
        d: Hamming distance threshold. Sequences will be corrected to the whitelist
            with hamming distance <= ``d``. Defaults to 1.
        n_threads: Number of threads to use. Defaults to 1.
        show_progress: Whether to display a progress bar. Defaults to True.

    Raises:
        SequenceError: If all the lengths of each sequence, qualitiy and
            whitelisted sequence are not equal, the number of sequences and
            qualities provided are not equal or the whitelist contains duplicates.

    Returns:
        The corrections as a dict of sequence to correction mappings. Note that
        the return type is different from :func:`correct_sequences_to_whitelist`.
    """
    # Check number of sequences and their lengths match with provided qualities
    if len(set(whitelist)) != len(whitelist):
        raise SequenceError('`whitelist` contains duplicates')
    for seq in sequences:
        if len(seq) != len(sequences[0]):
            raise SequenceError(
                'all sequences in `sequences` must be of same length'
            )
    for bc in whitelist:
        if len(bc) != len(whitelist[0]):
            raise SequenceError(
                'all sequences in `whitelist` must be of same length'
            )
    if len(sequences[0]) != len(whitelist[0]):
        raise SequenceError(
            'all sequences in `sequences` and `whitelist` must be of same length'
        )

    corrections = {}

    # Step 1: find exact matches first
    sequences_set = set(sequences)
    whitelist_set = set(whitelist)
    unmatched_sequences = []
    for seq in progress(list(sequences_set), desc='[1/3] Finding exact matches',
                        disable=not show_progress):
        if seq in whitelist_set:
            corrections[seq] = seq
        else:
            unmatched_sequences.append(seq)

    # Step 2: Construct whitelist mask
    whitelist_arrays = np.zeros(
        (len(whitelist), len(NUCLEOTIDES_STRICT), len(whitelist[0])),
        dtype=bool
    )
    for i, bc in enumerate(progress(whitelist, desc='[2/3] Constructing masks',
                                    disable=not show_progress)):
        whitelist_arrays[i] = _sequence_to_array(bc)

    # Step 3: Find all mismatch masks for which hamming distance <= d
    for i, (indices, masks) in enumerate(utils.ParallelWithProgress(
            n_jobs=n_threads, total=len(unmatched_sequences),
            desc='[3/3] Finding mismatches', disable=not show_progress
    )(delayed(_mismatch_masks)(_sequence_to_array(seq), whitelist_arrays, d=d)
      for seq in unmatched_sequences)):
        indices = np.array(indices, dtype=int)
        masks = np.array(masks, dtype=bool)
        seq = unmatched_sequences[i]

        if indices.shape[0] == 0:
            corrections[seq] = None
            continue

        # At this point, all exact matches have been processed, so we don't
        # need to worry about those.
        match_indices = indices[masks.sum(axis=1) == 0]
        if len(match_indices) == 1:
            corrections[seq] = whitelist[match_indices[0]]
        elif indices.shape[0] == 1:
            corrections[seq] = whitelist[indices[0]]
        else:
            corrections[seq] = None

    return corrections
