import gzip
import os
from abc import abstractmethod
from typing import Optional

from joblib import Parallel
from tqdm import tqdm


class ParallelWithProgress(Parallel):
    """Wrapper around joblib.Parallel that uses tqdm to print execution progress.
    Taken from https://stackoverflow.com/a/61900501
    """

    def __init__(
        self,
        use_tqdm: bool = True,
        total: Optional[int] = None,
        desc: Optional[str] = None,
        *args,
        **kwargs
    ):
        self._use_tqdm = use_tqdm
        self._total = total
        self._desc = desc
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        with tqdm(disable=not self._use_tqdm, total=self._total, smoothing=0,
                  desc=self._desc) as self._pbar:
            return Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()


def is_gzip(path: str):
    if os.path.isfile(path):
        magic = b'\x1f\x8b'
        with open(path, 'rb') as f:
            return magic == f.read(len(magic))
    else:
        return path.endswith('.gz')


def open_as_text(path: str, mode: str):
    return gzip.open(path, f'{mode}t') if is_gzip(path) else open(path, mode)


def all_exists(*paths):
    return all(os.path.exists(path) for path in paths)


class FileWrapper:

    def __init__(self, path: str, mode: str = 'r'):
        self.path = path
        self.mode = mode
        self.fp = None
        self.closed = False

        # Immediately open file descriptor
        self._open()

    @property
    def is_gzip(self):
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
        self.fp = open_as_text(self.path, self.mode)

    def close(self):
        if not self.fp.closed:
            self.fp.close()

    def reset(self):
        self.close()
        self.__init__(self.path, self.mode)

    def tell(self):
        return self.fp.tell()

    @abstractmethod
    def read(self):
        pass

    @abstractmethod
    def write(self, entry):
        pass
