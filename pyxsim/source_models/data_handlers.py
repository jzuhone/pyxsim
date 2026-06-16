import numpy as np


class YTDataHandler:
    def return_chunk(self, source_model, chunk):
        self.out_chunk = source_model.return_chunk(chunk)

    def __getitem__(self, key):
        if key[1] not in self.out_chunk:
            self.out_chunk[key] = np.ravel(self.in_chunk[key].d)
        return self.out_chunk[key[1]]
