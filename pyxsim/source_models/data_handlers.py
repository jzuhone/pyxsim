import numpy as np


class YTDataHandler:
    def __init__(self, in_chunk):
        self.in_chunk = in_chunk
        self.out_chunk = {}

    def __getitem__(self, key):
        if key[1] not in self.out_chunk:
            self.out_chunk[key] = np.ravel(self.in_chunk[key].d)
        return self.out_chunk[key[1]]
