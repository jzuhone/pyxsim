import pyxsim
import yt
from yt.utilities.cosmology import Cosmology
import numpy as np
from pyxsim.lib.sky_functions import pixel_to_cel
from pyxsim.utils import parse_value, mylog

axis_wcs = [[1,2],[0,2],[0,1]]


def make_grid_source(fn, axis, width, center, redshift, area,
                     exp_time, source_model, sky_center, fov,
                     simput_prefix, depth=None, cosmology=None, dist=None,
                     absorb_model=None, nH=None, no_shifting=False,
                     sigma_pos=None, kernel="top_hat", overwrite=False,
                     prng=None):
    from pyxsim.lib.sky_functions import pixel_to_cel

    sky_center = np.array(sky_center)

    ds = yt.load(fn)

    if cosmology is None:
        if hasattr(ds, 'cosmology'):
            cosmo = ds.cosmology
        else:
            cosmo = Cosmology()
    else:
        cosmo = cosmology

    axis = ds.coordinates.axis_id.get(axis, axis)
    center = ds.coordinates.sanitize_center(center, axis)[0]
    width = ds.coordinates.sanitize_width(axis, width, depth)
    xwidth = width[0].to("code_length")
    ywidth = width[1].to("code_length")
    if len(width) == 3:
        depth = width[2].to("code_length")
    else:
        depth = max(xwidth, ywidth)
    fov = parse_value(fov, "arcmin")

    if dist is None:
        fov_width = ds.quan(cosmo.angular_scale(0.0, redshift)*fov)
        fov_width.convert_to_units("code_length")
        D_A = ds.quan(cosmo.angular_diameter_distance(0.0, redshift))
        D_A.convert_to_units('code_length')
    else:
        D_A = parse_value(dist, "Mpc", ds=ds).to("code_length")
        fov_width = fov.to_value("radian")*D_A

    mylog.info("Linear width of field of view is {:.2f} kpc.".format(fov_width.to_value('kpc')))

    nx = int(np.ceil(xwidth / fov_width))
    ny = int(np.ceil(ywidth / fov_width))

    mylog.info("Number of FOV to cover source is {:d} x {:d})".format(nx, ny))

    axisx, axisy = axis_wcs[axis]

    outfile = "{}_grid.txt".format(simput_prefix)
    f = open(outfile, "w")
    f.write("# {}_simput.fits\n".format(simput_prefix))

    k = 0

    for i in range(nx):
        for j in range(ny):
            box_center = center.copy()
            box_center[axisx] += (i-0.5*(nx-1))*fov_width
            box_center[axisy] += (j-0.5*(ny-1))*fov_width
            le = box_center - 0.5*fov_width
            re = box_center + 0.5*fov_width
            le[axis] = box_center[axis] - 0.5*depth
            re[axis] = box_center[axis] + 0.5*depth
            box = ds.box(le, re)
            photons = pyxsim.PhotonList.from_data_source(box, redshift,
                                                         area, exp_time,
                                                         source_model,
                                                         center=box_center,
                                                         dist=dist,
                                                         cosmology=cosmo)
            xsky = np.array([(box_center-center)[axisx]/D_A])
            ysky = np.array([(box_center-center)[axisy]/D_A])
            pixel_to_cel(xsky, ysky, sky_center)
            ra = xsky[0]
            dec = ysky[0]
            events = photons.project_photons("xyz"[axis], (ra, dec),
                                             absorb_model=absorb_model,
                                             nH=nH, no_shifting=no_shifting,
                                             sigma_pos=sigma_pos, kernel=kernel,
                                             prng=prng)
            del photons

            phlist_prefix = "{:s}_{:d}_{:d}".format(simput_prefix, i, j)
            if k == 0:
                append = False
            else:
                append = True
            events.write_simput_file(phlist_prefix, overwrite=overwrite,
                                     append=append, simput_prefix=simput_prefix)

            del events

            phlist = "{}_phlist.fits".format(phlist_prefix)
            f.write("{:d}\t{:s}\t{:.2f}\t{:.2f}".format(k, phlist, ra, dec))
            k += 1

            f.flush()
            f.close()

            return outfile
