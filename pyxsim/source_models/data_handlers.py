from pyxsim.utils import DummyDataSource

try:
    from yt.data_objects.static_output import Dataset as YTDataset
except ImportError:
    YTDataset = DummyDataSource


class DataHandler:
    pass


class YTDataHandler(DataHandler):
    def __init__(self, data_source):
        self.data_source = data_source
        self.ds = data_source.ds

    def get_field_info(self, field):
        fi = self.ds._get_field_info(field)
        return fi

    def process_array(self, array, to_unit=None):
        if to_unit:
            return array.to_value(to_unit)
        else:
            return array.d


def find_data_handler(data_source):
    if isinstance(getattr(data_source, "ds", None), YTDataset):
        return YTDataHandler(data_source)
    else:
        raise NotImplementedError("Data source not supported")
