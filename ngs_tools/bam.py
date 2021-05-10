from typing import Any, Callable, Optional

import pysam
from tqdm import tqdm

from .fastq import Fastq, Read


class BamError(Exception):
    pass


def map_bam(
    bam_path: str,
    map_func: Callable[[pysam.AlignedSegment], Any],
    n_threads: int = 1
):
    """Generator to map an arbitrary function to every read and return its return
    values.

    Args:
        bam_path: Path to the BAM file
        map_func: Function that takes a :class:`pysam.AlignedSegment` object and
            returns some value
        n_threads: Number of threads to use. Defaults to 1.

    Yields:
        ``map_func`` applied to each read in the BAM file
    """
    with pysam.AlignmentFile(bam_path, 'rb', threads=n_threads,
                             check_sq=False) as f:
        for read in tqdm(f.fetch(until_eof=True), desc='Mapping BAM',
                         smoothing=0):
            yield map_func(read)


def apply_bam(
    bam_path: str,
    apply_func: Callable[[pysam.AlignedSegment],
                         Optional[pysam.AlignedSegment]],
    out_path: str,
    n_threads: int = 1,
):
    """Apply an arbitrary function to every read in a BAM. Reads for which the
    function returns `None` are not written to the output BAM.

    Args:
        bam_path: Path to the BAM file
        apply_func: Function that takes a :class:`pysam.AlignedSegment` object and
            optionally returns :class:`pysam.AlignedSegment` objects
        out_path: Path to output BAM file
        n_threads: Number of threads to use. Defaults to 1.

    Returns:
        Path to written BAM
    """
    with pysam.AlignmentFile(bam_path, 'rb', threads=n_threads,
                             check_sq=False) as f_in:
        with pysam.AlignmentFile(out_path, 'wb', template=f_in,
                                 threads=n_threads) as f_out:
            for read in tqdm(f_in.fetch(until_eof=True), desc='Applying BAM',
                             smoothing=0):
                result = apply_func(read)
                if result is not None:
                    f_out.write(result)
    return out_path


def tag_bam_with_fastq(
    bam_path: str,
    fastq_path: str,
    tag_func: Callable[[Read], dict],
    out_path: str,
    check_name: bool = True,
    n_threads: int = 1,
):
    """Add tags to BAM entries using sequences from a FASTQ file.

    Internally, this function calls :func:`apply_bam`.

    Args:
        bam_path: Path to the BAM file
        fastq_path: Path to FASTQ file
        tag_func: Function that takes a :class:`pysam.AlignedSegment` object and
            returns a dictionary of tags
        out_path: Path to output BAM file
        check_name: Whether or not to raise a :class:`BamError` if the FASTQ does not
            contain a read in the BAM
        n_threads: Number of threads to use. Defaults to 1.

    Returns:
        Path to written BAM
    """
    tags = {
        read.name: tag_func(read)
        for read in
        tqdm(Fastq(fastq_path), smoothing=0, desc='Extracting tags')
    }

    def apply_func(al: pysam.AlignedSegment):
        if al.query_name in tags:
            al.set_tags(list(tags[al.query_name].items()))
        elif check_name:
            raise BamError(f'Missing read `{al.query_name}` in FASTQ')
        return al

    return apply_bam(bam_path, apply_func, out_path, n_threads)


def filter_bam(
    bam_path: str,
    filter_func: Callable[[pysam.AlignedSegment], bool],
    out_path: str,
    n_threads: int = 1,
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

    Returns:
        Path to written BAM
    """
    return apply_bam(
        bam_path, lambda al: None
        if filter_func(al) is False else al, out_path, n_threads
    )
