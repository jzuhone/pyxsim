import numpy as np
from yt.units.yt_array import YTQuantity
from six import string_types
from collections import defaultdict
import h5py
import os
from astropy.units import Quantity
import logging
from more_itertools import always_iterable

try:
    # yt 3.x
    from yt.funcs import issue_deprecation_warning
except ImportError:
    # yt 4.x
    from yt._maintenance.deprecation import issue_deprecation_warning

pyxsimLogger = logging.getLogger("pyxsim")

ufstring = "%(name)-3s : [%(levelname)-9s] %(asctime)s %(message)s"
cfstring = "%(name)-3s : [%(levelname)-18s] %(asctime)s %(message)s"

pyxsim_sh = logging.StreamHandler()
# create formatter and add it to the handlers
formatter = logging.Formatter(ufstring)
pyxsim_sh.setFormatter(formatter)
# add the handler to the logger
pyxsimLogger.addHandler(pyxsim_sh)
pyxsimLogger.setLevel('INFO')
pyxsimLogger.propagate = False

mylog = pyxsimLogger


def force_unicode(value):
    if hasattr(value, 'decode'):
        return value.decode('utf8')
    else:
        return value


def parse_value(value, default_units, ds=None):
    if isinstance(value, Quantity):
        value = YTQuantity.from_astropy(value)
    if ds is None:
        quan = YTQuantity
    else:
        quan = ds.quan
    if isinstance(value, YTQuantity):
        return quan(value.v, value.units).in_units(default_units)
    elif isinstance(value, tuple):
        return quan(value[0], value[1]).in_units(default_units)
    else:
        return quan(value, default_units)


def validate_parameters(first, second, skip=[]):
    from h5py import Dataset
    keys1 = list(first.keys())
    keys2 = list(second.keys())
    keys1.sort()
    keys2.sort()
    if keys1 != keys2:
        raise RuntimeError("The two inputs do not have the same parameters!")
    for k1, k2 in zip(keys1, keys2):
        if k1 not in skip:
            if isinstance(first[k1], Dataset):
                v1 = first[k1][()]
            else:
                v1 = first[k1]
            if isinstance(first[k2], Dataset):
                v2 = first[k2][()]
            else:
                v2 = first[k2]
            if isinstance(v1, (string_types, bytes)) or isinstance(v2, (string_types, bytes)):
                check_equal = v1 == v2
            else:
                check_equal = np.allclose(np.array(v1), np.array(v2), rtol=0.0, atol=1.0e-10)
            if not check_equal:
                raise RuntimeError("The values for the parameter '%s' in the two inputs" % k1 +
                                   " are not identical (%s vs. %s)!" % (v1, v2))


def merge_files(input_files, output_file, overwrite=False,
                add_exposure_times=False):
    r"""
    Helper function for merging PhotonList or EventList HDF5 files.

    Parameters
    ----------
    input_files : list of strings
        List of filenames that will be merged together.
    output_file : string
        Name of the merged file to be outputted.
    overwrite : boolean, default False
        If a the output file already exists, set this to True to
        overwrite it.
    add_exposure_times : boolean, default False
        If set to True, exposure times will be added together. Otherwise,
        the exposure times of all of the files must be the same.

    Examples
    --------
    >>> from pyxsim import merge_files
    >>> merge_files(["events_0.h5","events_1.h5","events_3.h5"], "events.h5",
    ...             overwrite=True, add_exposure_times=True)

    Notes
    -----
    Currently, to merge files it is mandated that all of the parameters have the
    same values, with the exception of the exposure time parameter "exp_time". If
    add_exposure_times=False, the maximum exposure time will be used.
    """
    if os.path.exists(output_file) and not overwrite:
        raise IOError("Cannot overwrite existing file %s. " % output_file +
                      "If you want to do this, set overwrite=True.")

    f_in = h5py.File(input_files[0], "r")
    f_out = h5py.File(output_file, "w")

    exp_time_key = ""
    p_out = f_out.create_group("parameters")
    for key, param in f_in["parameters"].items():
        if key.endswith("exp_time"):
            exp_time_key = key
        else:
            p_out[key] = param[()]

    skip = [exp_time_key] if add_exposure_times else []
    for fn in input_files[1:]:
        f = h5py.File(fn, "r")
        validate_parameters(f_in["parameters"], f["parameters"], skip=skip)
        f.close()

    f_in.close()

    data = defaultdict(list)
    tot_exp_time = 0.0

    for i, fn in enumerate(input_files):
        f = h5py.File(fn, "r")
        if add_exposure_times:
            tot_exp_time += f["/parameters"][exp_time_key][()]
        else:
            tot_exp_time = max(tot_exp_time, f["/parameters"][exp_time_key][()])
        for key in f["/data"]:
            data[key].append(f["/data"][key][:])
        f.close()

    p_out[exp_time_key] = tot_exp_time

    d = f_out.create_group("data")
    for k in data:
        d.create_dataset(k, data=np.concatenate(data[k]))

    f_out.close()

key_warning = "As of pyXSIM v2.0.0 this key to the %s has been deprecated, " + \
              "and will be removed in a future release. Use the '%s' key instead."


class ParameterDict(object):
    def __init__(self, param_dict, list_type, old_keys):
        self.param_dict = param_dict
        self.list_type = list_type
        self.old_keys = old_keys

    def __setitem__(self, key, value):
        self.param_dict[key] = value

    def __getitem__(self, key):
        if key in self.old_keys:
            k = self.old_keys[key]
            mylog.warning(key_warning % ("%s parameters" % self.list_type, k))
        else:
            k = key
        return self.param_dict[k]

    def keys(self):
        return self.param_dict.keys()

    def values(self):
        return self.param_dict.values()

    def items(self):
        return self.param_dict.items()

    def __contains__(self, key):
        if key in self.old_keys:
            mylog.warning(key_warning % ("%s parameters" % self.list_type, self.old_keys[key]))
            return True
        return key in self.param_dict

    def __repr__(self):
        return self.param_dict.__repr__()


def ensure_list(obj):
    return list(always_iterable(obj))