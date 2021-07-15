from functools import partial

from tqdm import tqdm

progress = partial(tqdm, ascii=True, smoothing=0.1)
