import h5py
import numpy as np
from tqdm.auto import tqdm
from yt.visualization.volume_rendering.off_axis_projection import off_axis_projection

from pyxsim.utils import get_normal_and_north, parse_value


def make_absorption_map(
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

    L, north_vector, orient = get_normal_and_north(normal, north_vector=north_vector)

    width = parse_value(width, "kpc", ds)
    depth = parse_value(depth, "kpc", ds)

    nH = np.zeros((nwidth, nwidth, ndepth))

    pbar = tqdm(
        leave=True,
        total=ndepth,
        desc="Determining a cube of neutral hydrogen column density ",
    )

    w = ds.coordinates.sanitize_width(normal, width, depth)
    dz = w[2] / ndepth

    if isinstance(normal, str):
        id = ds.coordinates.axis_id[normal]
        le = center - 0.5 * w
        re = center + 0.5 * w
        for i in range(ndepth):
            le[id] += (i + 1) * dz
            box = ds.box(le, re)
            prj = ds.proj(field, normal, center=center, data_source=box)
            frb = prj.to_frb(width, nwidth)
            nH[:, :, i] = frb[field].d
            pbar.update()
    else:
        bw = ds.arr(w)
        re = center + 0.5 * depth * L
        for i in range(ndepth):
            bw[2] = (ndepth - i) * dz
            bc = re - 0.5 * (ndepth - i) * dz * L
            dk = ds.disk(bc, L, bw[0], 0.5 * bw[2])
            img = off_axis_projection(
                dk, center, L, bw, (nwidth, nwidth), field, north_vector=north_vector
            )
            nH[:, :, i] = np.asarray(img)
            pbar.update()

    pbar.close()

    wbins = np.linspace(-0.5 * width.d, 0.5 * width.d, nwidth + 1)
    dbins = np.linspace(-0.5 * depth.d, 0.5 * depth.d, ndepth + 1)
    wmid = 0.5 * (wbins[1:] + wbins[:-1])
    dmid = 0.5 * (dbins[1:] + dbins[:-1])
    with h5py.File(outfile, "w") as f:
        p = f.create_group("parameters")
        p.attrs["normal"] = L
        p.attrs["north"] = north_vector
        d = f.create_group("data")
        d.create_dataset("wmid", data=wmid)
        d.create_dataset("dmid", data=dmid)
        d.create_dataset("nH", data=nH * 1.0e-22)
