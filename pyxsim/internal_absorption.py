import h5py
import numpy as np
from tqdm.auto import tqdm
from yt.config import ytcfg
from yt.utilities.logger import ytLogger
from yt.visualization.volume_rendering.off_axis_projection import off_axis_projection

from pyxsim.utils import get_normal_and_north, parse_value


def make_column_density_map(
    ds,
    normal,
    center,
    width,
    depth,
    nwidth,
    ndepth,
    outfile,
    field=("gas", "H_p0_number_density"),
    north_vector=None,
):
    """
    Create a cube of neutral hydrogen column density to be used in the
    absorption of photons.

    Parameters
    ----------
    ds : yt Dataset
        The dataset to use when creating the column density cube.
    normal : character or array-like
        Normal vector to the plane of projection. If "x", "y", or "z", will
        assume to be along that axis (and will probably be faster). Otherwise,
        should be an off-axis normal vector, e.g [1.0, 2.0, -3.0]
    center : string or array_like
        The origin of the photon spatial coordinates. Accepts "c", "max", or
        a coordinate. If array-like and without units, it is assumed to be in
        units of kpc.
    width : float, tuple, or unyt_quantity
        The width of the cube in kpc. If a float, will assume units of kpc.
    depth : float, tuple, or unyt_quantity
        The width of the cube in kpc. If a float, will assume units of kpc.
    nwidth : integer
        The number of cells on a side in both of the sky directions of the
        hydrogen column density cube.
    ndepth : integer
        The number of cells on a side along the depth of the column density
        cube.
    outfile : string
        The HDF5 file to be written containing the cube.
    field : 2-tuple of strings, optional
        The yt field to be used for the neutral hydrogen density. Default is
        ("gas", "H_p0_number_density").
    north_vector : a sequence of floats
        A vector defining the "up" direction. This option sets the
        orientation of the plane of projection. If not set, an arbitrary
        grid-aligned north_vector perpendicular to the normal is chosen.
        Ignored in the case where a particular axis (e.g., "x", "y", or
        "z") is explicitly specified.
    """
    L, north_vector, orient = get_normal_and_north(normal, north_vector=north_vector)

    width = parse_value(width, "kpc", ds)
    depth = parse_value(depth, "kpc", ds)

    nH = np.zeros((nwidth, nwidth, ndepth))

    pbar = tqdm(
        leave=True,
        total=ndepth,
        desc="Determining a cube of neutral hydrogen column density ",
    )

    w = ds.arr(ds.coordinates.sanitize_width(normal, width, depth))
    center, _ = ds.coordinates.sanitize_center(center, normal)
    dz = w[2] / ndepth

    if isinstance(normal, str):
        le = center.copy()
        re = center.copy()
        dirs = [
            ds.coordinates.x_axis[normal],
            ds.coordinates.y_axis[normal],
            ds.coordinates.axis_id[normal],
        ]
        for ii, ax in enumerate(dirs):
            le[ax] = center[ax] - 0.5 * w[ii]
            re[ax] = center[ax] + 0.5 * w[ii]
        id = dirs[2]
        lei = le.copy()
        old_level = int(ytcfg.get("yt", "log_level"))
        ytLogger.setLevel(40)
        for i in range(ndepth):
            lei[id] = le[id] + i * dz
            box = ds.box(lei, re)
            prj = ds.proj(field, normal, center=center, data_source=box)
            frb = prj.to_frb(width, nwidth)
            nH[:, :, i] = frb[field].d
            pbar.update()
        ytLogger.setLevel(old_level)
    else:
        re = center + 0.5 * depth * L
        for i in range(ndepth):
            w[2] = (ndepth - i) * dz
            bc = re - 0.5 * (ndepth - i) * dz * L
            dk = ds.disk(bc, L, w[0], 0.5 * w[2])
            img = off_axis_projection(
                dk, center, L, w, (nwidth, nwidth), field, north_vector=north_vector
            )
            nH[:, :, i] = np.asarray(img)
            pbar.update()

    pbar.close()

    wbins = np.linspace(-0.5 * width.d, 0.5 * width.d, nwidth + 1)
    dbins = np.linspace(-0.5 * depth.d, 0.5 * depth.d, ndepth + 1)
    with h5py.File(outfile, "w") as f:
        p = f.create_group("parameters")
        p.attrs["normal"] = L
        p.attrs["north"] = north_vector
        d = f.create_group("data")
        d.create_dataset("wbins", data=wbins)
        d.create_dataset("dbins", data=dbins)
        d.create_dataset("nH", data=nH * 1.0e-22)
