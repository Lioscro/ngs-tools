from collections import Counter
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import pysam

from .fastq import Fastq, Read
from .progress import progress


class BamError(Exception):
    pass


def map_bam(
    bam_path: str,
    map_func: Callable[[pysam.AlignedSegment], Any],
    n_threads: int = 1,
    show_progress: bool = False,
):
    """Generator to map an arbitrary function to every read and return its return
    values.

    Args:
        bam_path: Path to the BAM file
        map_func: Function that takes a :class:`pysam.AlignedSegment` object and
            returns some value
        n_threads: Number of threads to use. Defaults to 1.
        show_progress: Whether to display a progress bar. Defaults to False.

    Yields:
        ``map_func`` applied to each read in the BAM file
    """
    with pysam.AlignmentFile(bam_path, 'rb', threads=n_threads,
                             check_sq=False) as f:
        for read in progress(f.fetch(until_eof=True), desc='Mapping BAM',
                             disable=not show_progress):
            yield map_func(read)


def apply_bam(
    bam_path: str,
    apply_func: Callable[[pysam.AlignedSegment],
                         Optional[pysam.AlignedSegment]],
    out_path: str,
    n_threads: int = 1,
    show_progress: bool = False,
):
    """Apply an arbitrary function to every read in a BAM. Reads for which the
    function returns `None` are not written to the output BAM.

    Args:
        bam_path: Path to the BAM file
        apply_func: Function that takes a :class:`pysam.AlignedSegment` object and
            optionally returns :class:`pysam.AlignedSegment` objects
        out_path: Path to output BAM file
        n_threads: Number of threads to use. Defaults to 1.
        show_progress: Whether to display a progress bar. Defaults to False.

    Returns:
        Path to written BAM
    """
    with pysam.AlignmentFile(bam_path, 'rb', threads=n_threads,
                             check_sq=False) as f_in:
        with pysam.AlignmentFile(out_path, 'wb', template=f_in,
                                 threads=n_threads) as f_out:
            for read in progress(f_in.fetch(until_eof=True),
                                 desc='Applying BAM',
                                 disable=not show_progress):
                result = apply_func(read)
                if result is not None:
                    f_out.write(result)
    return out_path


def count_bam(
    bam_path: str,
    filter_func: Optional[Callable[[pysam.AlignedSegment], bool]] = None,
    n_threads: int = 1,
    show_progress: bool = False,
) -> int:
    """Count the number of BAM entries. Optionally, a function may be provided to
    only count certain alignments.

    Args:
        bam_path: Path to BAM
        filter_func: Function that takes a :class:`pysam.AlignedSegment` object and
            returns True for reads to be counted and False otherwise
        n_threads: Number of threads to use. Defaults to 1.
        show_progress: Whether to display a progress bar. Defaults to False.

    Returns:
        Number of alignments in BAM
    """
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False,
                             threads=n_threads) as f:
        try:
            n = 0
            for index in f.get_index_statistics():
                n += index.total
            return n
        except ValueError:
            pass
    return sum(
        map_bam(
            bam_path,
            filter_func or (lambda al: True),
            n_threads=n_threads,
            show_progress=show_progress
        )
    )


def split_bam(
    bam_path: str,
    split_prefix: str,
    split_func: Optional[Callable[[pysam.AlignedSegment], str]] = None,
    n: Optional[int] = None,
    n_threads: int = 1,
    check_pair_groups: bool = True,
    show_progress: bool = False,
) -> Dict[str, Tuple[str, int]]:
    """Split a BAM into many parts, either by the number of reads or by an
    arbitrary function. Only one of ``split_func`` or ``n`` must be provided.
    Read pairs are always written to the same file.

    This function makes two passes through the BAM file. The first pass is to
    identify which reads must be written together (i.e. are pairs). The second
    pass is to actually extract the reads and write them to the appropriate split.

    The following procedure is used to identify pairs.
    1) The ``.is_paired`` property is checked to be True.
    2) If the read is uanligned, at most one other unaligned read with the same
       read name is allowed to be in the BAM. This other read is its mate.
       If the read is aligned, it should have the ``HI`` BAM tag indicating the
       alignment index. If no ``HI`` tag is present, then it is assumed only
       one alignment should be present for each read pair. If any of these
       constraints are not met, an exception is raised.

    Args:
        bam_path: Path to the BAM file
        split_prefix: File path prefix to all the split BAMs
        split_func: Function that takes a :class:`pysam.AlignedSegment` object and
            returns a string ID that is used to group reads into splits. All reads
            with a given ID will be written to a single BAM. Defaults to None.
        n: Number of BAMs to split into. Defaults to None.
        n_threads: Number of threads to use. Only affects reading. Writing is
            still serialized. Defaults to 1.
        check_pair_groups: When using ``split_func``, make sure that paired reads
            are assigned the same ID (and thus are split into the same BAM).
            Defaults to True.
        show_progress: Whether to display a progress bar. Defaults to False.

    Returns:
        Dictionary of tuples, where the first element is the path to a split BAM,
        and the second element is the number of BAM entries written to that split.
        The keys are either the string ID of each split (if ``split_func`` is used)
        or the split index (if ``n`` is used), and the values are paths.

    Raises:
        BamError: If any pair constraints are not met.
    """
    if (split_func is None) == (n is None):
        raise BamError('Exactly one of `split_func` or `n` must be provided.')

    n_reads = 0
    read_indices = []

    # Temporary dict that hold paired reads that we have only found one mate for
    pairs_cache = {}
    # Group ids for every read in order they appear in the BAM
    group_ids = []
    # Pairs that have been added
    pairs_added = set()

    # First, determine which reads should go to which split.
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False,
                             threads=n_threads) as f:
        for i, read in enumerate(f.fetch(until_eof=True)):
            n_reads += 1
            read_name = read.query_name

            group_id = None if not split_func else split_func(read)
            group_ids.append(group_id)

            if read.is_paired:
                alignment_index = None if read.is_unmapped or not read.has_tag(
                    'HI'
                ) else read.get_tag('HI')
                key = (read_name, alignment_index)
                if key in pairs_added:
                    raise BamError(
                        f'Found alignment for paired read {read_name} multiple times with '
                        f'index {alignment_index} (`None` means unmapped).'
                    )

                if key in pairs_cache:
                    mate_i = pairs_cache[key]
                    if split_func and check_pair_groups and group_id != group_ids[
                            mate_i]:
                        raise BamError(
                            f'Read {read_name} pairs are assigned different split IDs. '
                            'Use `check_pair_groups=False` to turn off this check.'
                        )
                    read_indices.append((mate_i, i))
                    del pairs_cache[key]
                    pairs_added.add(key)
                else:
                    pairs_cache[key] = i
            else:
                read_indices.append(i)

    # Make sure all paired reads have been covered
    if pairs_cache:
        raise BamError('Some paired reads do not have mates.')

    # When just splitting by using n, group_id is just an incrementing integer
    if not split_func:
        n_split = n_reads // n
        n_current = 0
        current_group = 0
        for idx in read_indices:
            if isinstance(idx, tuple):
                group_ids[idx[0]] = group_ids[idx[1]] = str(current_group)
                n_current += 2
            else:
                group_ids[idx] = str(current_group)
                n_current += 1

            if n_current > n_split:
                n_current = 0
                current_group += 1

    # Open a new file for each split.
    split_paths = {
        group_id: f'{split_prefix}_{group_id}.bam'
        for group_id in set(group_ids)
    }
    split_counts = Counter(group_ids)
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False,
                             threads=n_threads) as f:
        split_outs = {}
        try:
            for group_id, path in split_paths.items():
                split_outs[group_id] = pysam.AlignmentFile(
                    path, 'wb', template=f, threads=n_threads
                )

            for i, read in progress(enumerate(f.fetch(until_eof=True)),
                                    total=n_reads, desc='Splitting BAM',
                                    disable=not show_progress):
                split_outs[group_ids[i]].write(read)
        finally:
            for out in split_outs.values():
                out.close()
    return {
        group_id: (path, split_counts[group_id])
        for group_id, path in split_paths.items()
    }


def tag_bam_with_fastq(
    bam_path: str,
    fastq_path: Union[str, List[str]],
    tag_func: Union[Callable[[Read], dict], List[Callable[[Read], dict]]],
    out_path: str,
    check_name: bool = True,
    n_threads: int = 1,
    show_progress: bool = False,
):
    """Add tags to BAM entries using sequences from one or more FASTQ files.

    Internally, this function calls :func:`apply_bam`.

    Note:
        The tag keys generated from `tag_func` must contain unique keys of at
        most 2 characters.

    Args:
        bam_path: Path to the BAM file
        fastq_path: Path to FASTQ file. This option may be a list to extract
            tags from multiple FASTQ files. In this case, `tag_func` must also
            be a list of functions.
        tag_func: Function that takes a :class:`ngs_tools.fastq.Read` object and
            returns a dictionary of tags. When multiple FASTQs are being parsed
            simultaneously, each function needs to produce unique dictionary keys.
            Additionally, BAM tag keys may only be at most 2 characters. However,
            neither of these conditions are checked in favor of runtime.
        out_path: Path to output BAM file
        check_name: Whether or not to raise a :class:`BamError` if the FASTQ does not
            contain a read in the BAM
        n_threads: Number of threads to use. Defaults to 1.
        show_progress: Whether to display a progress bar. Defaults to False.

    Returns:
        Path to written BAM

    Raises:
        BamError: If only one of `fastq_path` and `tag_func` is a list, if
            both are lists but they have different lengths, if `check_name=True`
            but there are missing tags.
    """
    # Check arguments.
    if isinstance(fastq_path, list) and isinstance(tag_func, list):
        if len(fastq_path) != len(tag_func):
            raise BamError(
                '`fastq_path` and `tag_func` must contain the same number of elements.'
            )
        fastq_paths = fastq_path
        tag_funcs = tag_func
    elif not isinstance(fastq_path, list) and not isinstance(tag_func, list):
        fastq_paths = [fastq_path]
        tag_funcs = [tag_func]
    else:
        raise BamError(
            'Both `fastq_path` and `tag_func` must be lists, or neither must '
            'be lists.'
        )

    tags = {}
    for path, func in zip(fastq_paths, tag_funcs):
        for read in progress(Fastq(path), desc=f'Extracting tags from {path}',
                             disable=not show_progress):
            tags.setdefault(read.name, {}).update(func(read))

    def apply_func(al: pysam.AlignedSegment):
        if al.query_name in tags:
            al.set_tags(list(tags[al.query_name].items()))
        elif check_name:
            raise BamError(f'Missing read `{al.query_name}` in FASTQ')
        return al

    return apply_bam(
        bam_path, apply_func, out_path, n_threads, show_progress=show_progress
    )


def filter_bam(
    bam_path: str,
    filter_func: Callable[[pysam.AlignedSegment], bool],
    out_path: str,
    n_threads: int = 1,
    show_progress: bool = False,
):
    """Filter a BAM by applying the given function to each :class:`pysam.AlignedSegment`
    object. When the function returns False, the read is not written to the output
    BAM.

    Internally, this function calls :func:`apply_bam`.

    Args:
        bam_path: Path to the BAM file
        filter_func: Function that takes a :class:`pysam.AlignedSegment` object and
            returns False for reads to be filtered out
        out_path: Path to output BAM file
        n_threads: Number of threads to use. Defaults to 1.
        show_progress: Whether to display a progress bar. Defaults to False.

    Returns:
        Path to written BAM
    """
    return apply_bam(
        bam_path,
        lambda al: None if filter_func(al) is False else al,
        out_path,
        n_threads,
        show_progress=show_progress
    )
