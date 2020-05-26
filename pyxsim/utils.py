from yt.funcs import iterable
from yt.units.yt_array import YTQuantity, YTArray
from unyt import unyt_array
import sys
from yt.config import ytcfg
from astropy.units import Quantity
import logging

pyxsimLogger = logging.getLogger("pyxsim")

if ytcfg.getboolean("yt", "stdoutStreamLogging"):
    stream = sys.stdout
else:
    stream = sys.stderr

level = min(max(ytcfg.getint("yt", "loglevel"), 0), 50)
ufstring = "%(name)-3s : [%(levelname)-9s] %(asctime)s %(message)s"
cfstring = "%(name)-3s : [%(levelname)-18s] %(asctime)s %(message)s"

pyxsim_sh = logging.StreamHandler(stream=stream)
# create formatter and add it to the handlers
formatter = logging.Formatter(ufstring)
pyxsim_sh.setFormatter(formatter)
# add the handler to the logger
pyxsimLogger.addHandler(pyxsim_sh)
pyxsimLogger.setLevel(level)
pyxsimLogger.propagate = False

mylog = pyxsimLogger


def parse_value(value, default_units, ds=None):
    if isinstance(value, Quantity):
        value = YTQuantity.from_astropy(value)
    if ds is None:
        quan = YTQuantity
    else:
        quan = ds.quan
    if isinstance(value, YTQuantity):
        return quan(value.v, value.units).in_units(default_units)
    elif iterable(value):
        return quan(value[0], value[1]).in_units(default_units)
    else:
        return quan(value, default_units)


def isunitful(a):
    if isinstance(a, (YTArray, Quantity, unyt_array)):
        return True
    elif isinstance(a, tuple):
        try:
            unyt_array(a[0], a[1])
            return True
        except:
            pass
    return False
