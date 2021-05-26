import gzip
import os
import tempfile
from abc import abstractmethod
from typing import Any, Optional, TextIO

from joblib import Parallel
from tqdm import tqdm
from typing_extensions import Literal


class ParallelWithProgress(Parallel):
    """Wrapper around joblib.Parallel that uses tqdm to print execution progress.
    Taken from https://stackoverflow.com/a/61900501
    """

    def __init__(
        self,
        pbar: Optional[tqdm] = None,
        total: Optional[int] = None,
        desc: Optional[str] = None,
        *args,
        **kwargs
    ):
        self._pbar = pbar or tqdm(total=total, desc=desc, smoothing=0)
        super(ParallelWithProgress, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        try:
            return Parallel.__call__(self, *args, **kwargs)
        finally:
            self._pbar.close()

    def print_progress(self):
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()


def is_gzip(path: str) -> bool:
    """Check if a file is Gzipped by checking the magic string.

    Args:
        path: path to file

    Returns:
        True or False
    """
    if os.path.isfile(path):
        magic = b'\x1f\x8b'
        with open(path, 'rb') as f:
            return magic == f.read(len(magic))
    else:
        return path.endswith('.gz')


def open_as_text(path: str, mode: Literal['r', 'w']) -> TextIO:
    """Open a (possibly gzipped) file in text mode.

    Args:
        path: Path to file
        mode: Mode to open file in. Either ``r`` for read or ``w`` for write.

    Returns:
        Opened file pointer that supports ``read`` and ``write`` functions.
    """
    return gzip.open(path, f'{mode}t') if is_gzip(path) else open(path, mode)


def all_exists(*paths: str) -> bool:
    """Check whether all provided paths exist.

    Args:
        *paths: paths to files

    Returns:
        True if all files exist, False otherwise
    """
    return all(os.path.exists(path) for path in paths)


class FileWrapper:
    """Generic wrapper class for file-formats. Used to wrap file-format-specific
    implementations of reading and writing entries. This class is not designed to
    be initialized directly. Instead, it should be inherited by children that
    implements the ``read`` and ``write`` methods appropriately.

    The file is opened immediately as soon as the class is initialized. This class
    can also be used as a context manager to safely close the file pointer with a
    ``with`` block.

    Attributes:
        path: Path to the file
        mode: Open mode. Either ``r`` or ``w``.
        fp: File pointer
        closed: Whether the file has been closed
    """

    def __init__(self, path: str, mode: Literal['r', 'w'] = 'r'):
        """
        Args:
            path: Path to the file
            mode: Open mode. Either ``r`` or ``w``.
        """
        self.path = path
        self.mode = mode
        self.fp = None
        self.closed = False

        # Immediately open file descriptor
        self._open()

    @property
    def is_gzip(self) -> bool:
        """Whether or not the file is gzipped"""
        return is_gzip(self.path)

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.close()

    def __iter__(self):
        while True:
            try:
                yield self.read()
            except StopIteration:
                return

    def _open(self):
        """Open the file"""
        self.fp = open_as_text(self.path, self.mode)

    def close(self):
        """Close the (possibly already-closed) file"""
        if not self.fp.closed:
            self.fp.close()

    def reset(self):
        """Reset this wrapper by first closing the file and re-running initialization,
        which re-opens the file."""
        self.close()
        self.__init__(self.path, self.mode)

    def tell(self) -> int:
        """Get the current location of the file pointer"""
        return self.fp.tell()

    @abstractmethod
    def read(self) -> Any:
        """Read a single entry. This method must be overridden by children."""
        pass

    @abstractmethod
    def write(self, entry: Any):
        """Write a single entry. This method must be overridden by children."""
        pass


def mkstemp(dir: Optional[str] = None, delete: bool = False):
    """Wrapper for :func:`tempfile.mkstemp` that automatically closes the OS-level
    file descriptor. This function behaves like :func:`tempfile.mkdtemp` but for
    files.

    Args:
        dir: Directory to create the temporary file. This value is passed as
            the ``dir`` kwarg of :func:`tempfile.mkstemp`. Defaults to None.
        delete: Whether to delete the temporary file before returning.
            Defaults to False.

    Returns:
        path to the temporary file
    """
    fd, path = tempfile.mkstemp(dir=dir)
    os.close(fd)
    if delete:
        os.remove(path)
    return path
