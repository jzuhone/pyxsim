from unyt import unyt_array, unyt_quantity
from astropy.units import Quantity
import logging
from more_itertools import always_iterable

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
