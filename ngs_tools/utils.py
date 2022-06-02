import functools
import gzip
import io
import os
import pickle
import queue
import shutil
import stat
import subprocess
import tempfile
import threading
import time
from abc import abstractmethod
from contextlib import contextmanager
from operator import add
from typing import (
    Any,
    Callable,
    Generator,
    Iterable,
    List,
    Optional,
    TextIO,
    Tuple,
    Union,
)
from urllib.parse import urlparse
from urllib.request import urlopen, urlretrieve

from joblib import Parallel
from tqdm import tqdm
from typing_extensions import Literal

from .logging import logger
from .progress import progress


class suppress_stdout_stderr:
    """A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).
    https://github.com/facebook/prophet/issues/223
    """

    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


def retry(
    function: Callable,
    retries: int,
    args: Optional[tuple] = None,
    kwargs: Optional[dict] = None,
    retry_every: Optional[int] = None,
    backoff: bool = False,
    exceptions: Optional[Tuple[Exception]] = None,
) -> Any:
    """Utility function to retry a function some number of times, with optional
    exponential backoff.

    Args:
        function: Function to retry
        retries: Number of times to retry
        args: Function arguments
        kwargs: Dictionary of keyword arguments
        retry_every: Time to wait in seconds between retries. Defaults to no wait time.
        backoff: Whether or not to exponential backoff between retries
        exceptions: Tuple of exceptions to expect. Defaults to all exceptions.

    Returns:
        Whatever ``function`` returns
    """
    args = args or tuple()
    kwargs = kwargs or {}
    exceptions = exceptions or Exception
    failed = 0
    while True:
        try:
            return function(*args, **kwargs)
        except exceptions:
            failed += 1
            if failed >= retries:
                raise
            if retry_every:
                time.sleep(retry_every)
                if backoff:
                    retry_every *= 2


def retry_decorator(
    retries: int,
    retry_every: Optional[int] = None,
    backoff: bool = False,
    exceptions: Optional[Tuple[Exception]] = None,
) -> Callable:
    """Function decorator to retry a function on exceptions.

    Args:
        retries: Number of times to retry
        retry_every: Time to wait in seconds between retries. Defaults to no wait time.
        backoff: Whether or not to exponential backoff between retries
        exceptions: Tuple of exceptions to expect. Defaults to all exceptions.
    """

    def decorator(func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return retry(
                func,
                retries,
                args,
                kwargs,
                retry_every=retry_every,
                backoff=backoff,
                exceptions=exceptions
            )

        return wrapper

    return decorator


def run_executable(
    command: List[str],
    stdin=None,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    wait: bool = True,
    stream: bool = True,
    quiet: bool = False,
    returncode: int = 0,
    alias: bool = True,
) -> Union[subprocess.Popen, Tuple[subprocess.Popen, List[str], List[str]]]:
    """Execute a single shell command.

    Args:
        command: A list representing a single shell command
        stdin: Object to pass into the ``stdin`` argument for :class:``subprocess.Popen``.
            Defaults to None
        stdout: Object to pass into the `stdout` argument for :class:``subprocess.Popen``
            Defaults to ``subprocess.PIPE``
        stderr: Object to pass into the `stderr` argument for :class:``subprocess.Popen``,
            Defaults to ``subprocess.PIPE``
        wait: Whether to wait until the command has finished. Defaults to True
        stream: Whether to stream the output to the command line. Defaults to True
        quiet: Whether to not display anything to the command line and not check the return code.
            Defaults to False
        returncode: The return code expected if the command runs as intended.
            Defaults to 0
        alias: Whether to use the basename of the first element of `command`.
            Defaults to `True`

    Returns:
        A tuple of (the spawned process, string printed to stdout,
        string printed to stderr) if `wait=True`. Otherwise, just the
        spawned process.

    Raises:
        subprocess.CalledProcessError: If not ``quiet`` and the process
            exited with an exit code != ``exitcode``
    """
    command = [str(c) for c in command]
    if not quiet:
        c = command.copy()
        if alias:
            c[0] = os.path.basename(c[0])
        logger.debug(' '.join(c))
    p = subprocess.Popen(
        command,
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
        universal_newlines=wait,
        bufsize=1 if wait else -1,
    )

    # Helper function to read from a pipe and put the output to a queue.
    def reader(pipe, qu, stop_event, name):
        while not stop_event.is_set():
            for _line in pipe:
                line = _line.strip()
                qu.put((name, line))

    # Wait if desired.
    if wait:
        stdout = ''
        stderr = ''
        out = []
        out_queue = queue.Queue()
        stop_event = threading.Event()
        stdout_reader = threading.Thread(
            target=reader,
            args=(p.stdout, out_queue, stop_event, 'stdout'),
            daemon=True
        )
        stderr_reader = threading.Thread(
            target=reader,
            args=(p.stderr, out_queue, stop_event, 'stderr'),
            daemon=True
        )
        stdout_reader.start()
        stderr_reader.start()

        while p.poll() is None:
            while not out_queue.empty():
                name, line = out_queue.get()
                if stream and not quiet:
                    logger.debug(line)
                out.append(line)
                if name == 'stdout':
                    stdout += f'{line}\n'
                elif name == 'stderr':
                    stderr += f'{line}\n'
            else:
                time.sleep(0.1)

        # Stop readers & flush queue
        stop_event.set()
        time.sleep(1)
        while not out_queue.empty():
            name, line = out_queue.get()
            if stream and not quiet:
                logger.debug(line)
            out.append(line)
            if name == 'stdout':
                stdout += f'{line}\n'
            elif name == 'stderr':
                stderr += f'{line}\n'

        if not quiet and p.returncode != returncode:
            logger.error('\n'.join(out))
            raise subprocess.CalledProcessError(p.returncode, ' '.join(command))

    return (p, stdout, stderr) if wait else p


class ParallelWithProgress(Parallel):
    """Wrapper around joblib.Parallel that uses tqdm to print execution progress.
    Taken from https://stackoverflow.com/a/61900501
    """

    def __init__(
        self,
        pbar: Optional[tqdm] = None,
        total: Optional[int] = None,
        desc: Optional[str] = None,
        disable: bool = False,
        *args,
        **kwargs
    ):
        self._pbar = pbar or progress(total=total, desc=desc, disable=disable)
        super(ParallelWithProgress, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        try:
            return Parallel.__call__(self, *args, **kwargs)
        finally:
            self._pbar.close()

    def print_progress(self):
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()


def is_remote(path: str) -> bool:
    """Check if a string is a remote URL.

    Args:
        path: string to check

    Returns:
        True or False
    """
    return bool(urlparse(path).scheme)


def is_gzip(path: str) -> bool:
    """Check if a file is Gzipped by checking the magic string.

    Args:
        path: path to file

    Returns:
        True or False
    """
    magic = b'\x1f\x8b'
    if is_remote(path):
        with urlopen(path) as f:
            return magic == f.read(len(magic))
    elif os.path.isfile(path):
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


def decompress_gzip(gzip_path: str, out_path: str) -> str:
    """Decompress a gzip file to provided file path.

    Args:
        gzip_path: Path to gzip file
        out_path: Path to decompressed file

    Returns:
        Path to decompressed file
    """
    with gzip.open(gzip_path, 'rb') as f, open(out_path, 'wb') as out:
        shutil.copyfileobj(f, out)
    return out_path


def compress_gzip(file_path: str, out_path: str) -> str:
    """Compress a file into gzip.

    Args:
        file_path: Path to file
        out_dir: Path to compressed file

    Returns:
        Path to compressed file
    """
    with open(file_path, 'rb') as f, gzip.open(out_path, 'wb') as out:
        shutil.copyfileobj(f, out)
    return out_path


def concatenate_files(*paths: str, out_path: str):
    """Concatenates an arbitrary number of files into one file.

    Args:
        *paths: An arbitrary number of paths to files
        out_path: Path to place concatenated file

    Returns:
        Path to concatenated file
    """
    with open(out_path, 'wb') as out:
        for path in paths:
            with open(path, 'rb') as f:
                shutil.copyfileobj(f, out)

    return out_path


def concatenate_files_as_text(*paths: str, out_path: str) -> str:
    """Concatenates an arbitrary number of files into one TEXT file.

    Only supports plaintext and gzip files.

    Args:
        *paths: An arbitrary number of paths to files
        out_path: Path to place concatenated file

    Returns:
        Path to concatenated file
    """
    with open(out_path, 'w') as out:
        for path in paths:
            with open_as_text(path, 'r') as f:
                for line in f:
                    if not line.isspace():
                        out.write(line.strip() + '\n')

    return out_path


class TqdmUpTo(tqdm):
    """Wrapper around :func:`tqdm` so that it can be used with :func:`urlretrieve`.
    https://github.com/tqdm/tqdm/blob/master/examples/tqdm_wget.py
    """

    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        return self.update(b * bsize - self.n)  # also sets self.n = b * bsize


def download_file(url: str, path: str) -> str:
    """Download a remote file to the provided path while displaying a progress bar.

    Args:
        url: Remote url
        path: Local path to download the file to

    Returns:
        Path to downloaded file
    """
    with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024,
                  miniters=1) as t:
        urlretrieve(url, filename=path, reporthook=t.update_to, data=None)
        t.total = t.n
    return path


@contextmanager
def stream_file(url: str, path: str) -> str:
    """A context manager that creates a FIFO file to use for piping remote files
    into processes. This function must be used as a context manager (the ``with``
    keyword) so that any exceptions in the streaming thread may be captured.

    This function spawns a new thread to download the remote file into a FIFO
    file object. FIFO file objects are only supported on unix systems.

    Args:
        url: Url to the file
        path: Path to place FIFO file

    Yields:
        Path to FIFO file

    Raises:
        OSError: If the operating system does not support FIFO

    """
    try:
        os.mkfifo(path)
    except AttributeError:
        raise OSError(
            'Operating system does not support streaming FIFO files. '
            'Download the file instead.'
        )
    logger.debug(f'Streaming {url} to {path}')
    t = threading.Thread(target=urlretrieve, args=(url, path), daemon=True)
    t.start()
    try:
        yield path
    finally:
        # There are two possible cases of getting here:
        # 1) The file has been completely streamed and processed
        # 2) An error occurred in the with block
        # In the case of 1), join should return immediately. In the case of 2),
        # we need a timeout because there isn't any point in fully downloading the
        # file.
        t.join(timeout=10)


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

        # Immediately open file descriptor
        self._open()

    @property
    def is_remote(self) -> bool:
        return is_remote(self.path)

    @property
    def is_gzip(self) -> bool:
        return is_gzip(self.path)

    @property
    def closed(self) -> bool:
        if self.fp is None:
            return True
        return self.fp.closed

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
        # Make sure to close the previous fp, if it is still open.
        self.close()

        if self.is_remote:
            if self.mode != 'r':
                raise IOError(
                    f'Remote URLs may not be opened in `{self.mode}` mode'
                )

            self.fp = urlopen(self.path)
            if self.is_gzip:
                self.fp = io.TextIOWrapper(
                    gzip.GzipFile(fileobj=self.fp, mode=self.mode)
                )
        else:
            self.fp = open_as_text(self.path, self.mode)

    def close(self):
        """Close the (possibly already-closed) file"""
        if not self.closed:
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


def write_pickle(obj: object, path: str, *args, **kwargs) -> str:
    """Pickle a Python object and compress with Gzip.

    Any additional arguments and keyword arguments are passed to :func:`pickle.dump`.

    Args:
        obj: Object to pickle
        path: Path to save pickle

    Returns:
        Saved pickle path
    """
    with gzip.open(path, 'wb') as f:
        pickle.dump(obj, f, *args, **kwargs)
    return path


def read_pickle(path: str) -> object:
    """Load a Python pickle that was compressed with Gzip.

    Args:
        path: Path to pickle

    Returns:
        Unpickled object
    """
    with gzip.open(path, 'rb') as f:
        return pickle.load(f)


def flatten_dictionary(
    d: dict,
    keys: Optional[tuple] = None
) -> Generator[Tuple[tuple, object], None, None]:
    """Generator that flattens the given dictionary into 2-element tuples
    containing keys and values. For nested dictionaries, the keys are
    appended into a tuple.

    Args:
        d: Dictionary to flatten
        keys: Previous keys, defaults to None. Used exclusively for recursion.

    Yields:
        Flattened dictionary as (keys, value)
    """
    keys = keys or tuple()
    for k, v in d.items():
        new_keys = keys + (k,)
        if isinstance(v, dict):
            yield from flatten_dictionary(v, new_keys)
        else:
            yield new_keys, v


def flatten_iter(it: Iterable) -> Generator[object, None, None]:
    """Generator that flattens the given iterable, except for strings.

    Args:
        lst: Iterable to flatten

    Yields:
        Flattened iterable elements
    """

    def is_iterable(i):
        if isinstance(i, str):
            return False
        try:
            iter(i)
            return True
        except TypeError:
            return False

    if not is_iterable(it):
        yield it

    for element in it:
        if not is_iterable(element):
            yield element
        else:
            yield from flatten_iter(element)


def merge_dictionaries(
    d1: dict,
    d2: dict,
    f: Callable[[object, object], object] = add,
    default: object = 0
) -> dict:
    """Merge two dictionaries, applying an arbitrary function `f` to duplicate keys.
    Dictionaries may be nested.

    Args:
        d1: First dictionary
        d2: Second dictionary
        f: Merge function. This function should take two arguments and return one,
            defaults to `+`
        default: Default value or callable to use for keys not present in either
            dictionary, defaults to `0`

    Returns:
        Merged dictionary
    """

    def setdefault_nested(d, t, value):
        inner = d
        for k in t[:-1]:
            inner = inner.setdefault(k, {})
        return inner.setdefault(t[-1], value)

    def get_nested(d, t, default=None):
        inner = d
        for k in t[:-1]:
            if k not in inner:
                return default() if callable(default) else default
            inner = inner[k]
        return inner.get(t[-1], default() if callable(default) else default)

    # Extract all keys
    d1_keys = [key for key, value in flatten_dictionary(d1)]
    d2_keys = [key for key, value in flatten_dictionary(d2)]
    keys = list(set(d1_keys + d2_keys))

    merged = {}
    for key in sorted(keys):
        setdefault_nested(
            merged, key,
            f(get_nested(d1, key, default), get_nested(d2, key, default))
        )

    return merged


def flatten_dict_values(d: dict) -> list:
    """Extract all values from a nested dictionary.

    Args:
        d: Nested dictionary from which to extract values from

    Returns:
        All values from the dictionary as a list
    """
    if isinstance(d, dict):
        flattened = []
        for k, v in d.items():
            if isinstance(v, dict):
                flattened.extend(flatten_dict_values(v))
            else:
                flattened.append(v)
        return flattened
    else:
        return [d]


def set_executable(path: str):
    """Set the permissions of a file to be executable.
    """
    st = os.stat(path)
    os.chmod(path, st.st_mode | stat.S_IEXEC)
