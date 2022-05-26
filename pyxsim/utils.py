from unyt import unyt_array, unyt_quantity
from astropy.units import Quantity
import logging
from more_itertools import always_iterable
import numpy as np
from soxs.constants import elem_names, \
    atomic_weights, abund_tables


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


def parse_value(value, default_units, ds=None):
    if isinstance(value, Quantity):
        value = unyt_quantity.from_astropy(value)
    if ds is None:
        quan = unyt_quantity
    else:
        quan = ds.quan
    if isinstance(value, unyt_quantity):
        return quan(value.v, value.units).in_units(default_units)
    elif isinstance(value, tuple):
        return quan(value[0], value[1]).in_units(default_units)
    else:
        return quan(value, default_units)


def isunitful(a):
    if isinstance(a, (Quantity, unyt_array)):
        return True
    elif isinstance(a, tuple):
        try:
            unyt_array(a[0], a[1])
            return True
        except:
            pass
    return False


def ensure_list(obj):
    return list(always_iterable(obj))


def validate_parameters(first, second, skip=None):
    if skip is None:
        skip = []
    keys1 = list(first.keys())
    keys2 = list(second.keys())
    keys1.sort()
    keys2.sort()
    if keys1 != keys2:
        raise RuntimeError("The two inputs do not have the same parameters!")
    for k1, k2 in zip(keys1, keys2):
        if k1 not in skip:
            v1 = first[k1][()]
            v2 = first[k2][()]
            if isinstance(v1, (str, bytes)) or isinstance(v2, (str, bytes)):
                check_equal = v1 == v2
            else:
                check_equal = np.allclose(np.array(v1), np.array(v2), rtol=0.0, atol=1.0e-10)
            if not check_equal:
                raise RuntimeError(f"The values for the parameter '{k1}' in the two inputs"
                                   f" are not identical ({v1} vs. {v2})!")


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
    from collections import defaultdict
    from pathlib import Path
    import h5py
    if Path(output_file).exists() and not overwrite:
        raise IOError(f"Cannot overwrite existing file {output_file}. "
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
        with h5py.File(fn, "r") as f:
            validate_parameters(f_in["parameters"], f["parameters"], skip=skip)

    f_in.close()

    data = defaultdict(list)
    tot_exp_time = 0.0

    for i, fn in enumerate(input_files):
        with h5py.File(fn, "r") as f:
            if add_exposure_times:
                tot_exp_time += f["/parameters"][exp_time_key][()]
            else:
                tot_exp_time = max(tot_exp_time, f["/parameters"][exp_time_key][()])
            for key in f["/data"]:
                data[key].append(f["/data"][key][:])

    p_out[exp_time_key] = tot_exp_time

    d = f_out.create_group("data")
    for k in data:
        d.create_dataset(k, data=np.concatenate(data[k]))

    f_out.close()


def compute_elem_mass_fraction(elem, abund_table="angr"):
    if isinstance(elem, str):
        elem = elem_names.index(elem)
    atable = abund_tables[abund_table]
    mZ = (atomic_weights[3:]*atable[3:]).sum()
    mE = atomic_weights[elem]*atable[elem]
    return mE/mZ


def compute_H_abund(abund_table):
    return atomic_weights[1]/(atomic_weights*abund_tables[abund_table]).sum()