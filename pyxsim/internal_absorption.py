import h5py
import numpy as np
from tqdm.auto import tqdm
from yt.visualization.volume_rendering import off_axis_projection

from pyxsim.utils import parse_value


def internal_absorption(
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

    if isinstance(normal, str):
        ax = "xyz".index(normal)
        L = np.zeros(3)
        L[ax] = 1.0
    else:
        L = np.array(normal)
    L /= np.sqrt(np.dot(L, L))

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
            re[id] = le[id] + (i + 1) * dz
            box = ds.box(le, re)
            prj = ds.proj(field, normal, center=box.center, data_source=box)
            frb = prj.to_frb(width, nwidth)
            nH[:, :, i] = frb[field].d
            pbar.update()
    else:
        bw = ds.quan(w)
        le = center - 0.5 * depth * L
        for i in range(ndepth):
            bw[i] = (i + 1) * dz
            bc = le + 0.5 * (i + 1) * dz * L
            img = off_axis_projection(
                ds, bc, L, bw, nwidth, field, north_vector=north_vector
            )
            nH[:, :, i] = np.asarray(img)
            pbar.update()

    pbar.close()

    wbins = np.linspace(-0.5 * width.d, 0.5 * width.d, nwidth + 1)
    dbins = np.linspace(-0.5 * depth.d, 0.5 * depth.d, ndepth + 1)
    with h5py.File(outfile, "w") as f:
        p = f.create_group("parameters")
        p.attrs["normal"] = L
        d = f.create_group("data")
        d.create_dataset("wbins", data=wbins)
        d.create_dataset("dbins", data=dbins)
        d.create_dataset("nH", data=nH * 1.0e-22)
