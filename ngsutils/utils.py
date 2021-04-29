from joblib import Parallel
from tqdm import tqdm


class ParallelWithProgress(Parallel):
    """Wrapper around joblib.Parallel that uses tqdm to print execution progress.
    Taken from https://stackoverflow.com/a/61900501
    """

    def __init__(self, use_tqdm=True, total=None, *args, **kwargs):
        self._use_tqdm = use_tqdm
        self._total = total
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        with tqdm(disable=not self._use_tqdm, total=self._total) as self._pbar:
            return Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        if self._total is None:
            self._pbar.total = self.n_dispatched_tasks
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()


def is_gzip(path):
    magic = b'\x1f\x8b'
    with open(path, 'rb') as f:
        return magic == f.read(len(magic))
