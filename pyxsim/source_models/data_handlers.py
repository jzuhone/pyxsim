from pyxsim.utils import DummyDataSource

try:
    from yt.data_objects.static_output import Dataset as YTDataset
except ImportError:
    YTDataset = DummyDataSource


class DataHandler:
    pass


class YTDataHandler(DataHandler):
    def __init__(self, ds):
        self.ds = ds

    def get_field_info(self, field):
        fi = self.ds._get_field_info(field)
        return fi

    def process_array(self, array, to_unit=None):
        if to_unit:
            return array.to_value(to_unit)
        else:
            return array.d


def find_data_handler(ds):
    if isinstance(ds, YTDataset):
        return YTDataHandler(ds)
    else:
        raise NotImplementedError("Data source not supported")
