from scipy.interpolate import InterpolatedUnivariateSpline
from yt.frontends.stream.api import load_uniform_grid
from yt.funcs import ensure_list
from yt.units.yt_array import YTArray
import numpy as np

def _create_dataset(radius, fields, r, ddims, width, bbox,
                    length_unit, mass_unit, time_unit, nprocs):

    if isinstance(radius, tuple):
        radius = YTArray(radius[0], radius[1]).to(length_unit).d
    elif isinstance(radius, YTArray):
        radius = radius.to(length_unit).d
    else:
        radius = YTArray(radius, length_unit).d

    data = {}
    for name, arr in fields.items():
        if isinstance(arr, tuple):
            arr = arr[0]
            units = arr[1]
        elif isinstance(arr, YTArray):
            arr = arr.d
            units = str(arr.units)
        else:
            units = "dimensionless"
        f = InterpolatedUnivariateSpline(radius, arr)
        data[name] = (f(r), units)

    bbox = np.array([[-0.5*width, 0.5*width]]*3)

    ds = load_uniform_grid(data, ddims, length_unit=length_unit,
                           time_unit=time_unit, mass_unit=mass_unit,
                           bbox=bbox, nprocs=nprocs)

    return ds

def create_spherical_dataset(radius, fields, ddims, width, length_unit="kpc",
                             time_unit="Gyr", mass_unit="Msun", nprocs=64):
    fields = ensure_list(fields)

    R = 0.5*width
    x, y, z = np.mgrid[-R:R:ddims[0]*1j,
                       -R:R:ddims[1]*1j,
                       -R:R:ddims[2]*1j]
    r = np.sqrt(x**2 + y**2 + z**2)

    bbox = np.array([[-R, R]]*3)

    return _create_dataset(radius, fields, r, ddims, width, bbox,
                           length_unit, mass_unit, time_unit, nprocs)

def create_cylindrical_dataset(radius, fields, ddims, width, height, 
                               length_unit="kpc", time_unit="Gyr", 
                               mass_unit="Msun", nprocs=64):
    fields = ensure_list(fields)

    R = 0.5*width
    H = 0.5*height
    x, y, z = np.mgrid[-R:R:ddims[0]*1j,
                       -R:R:ddims[1]*1j,
                       -H:H:ddims[2]*1j]
    r = np.sqrt(x**2 + y**2)

    bbox = np.array([[-R, R], [-R, R], [-H, H]])

    return _create_dataset(radius, fields, r, ddims, width, bbox,
                           length_unit, mass_unit, time_unit, nprocs)
