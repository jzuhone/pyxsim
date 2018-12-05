from yt.utilities.cosmology import Cosmology
import numpy as np
from pyxsim.lib.sky_functions import pixel_to_cel
from pyxsim.photon_list import PhotonList
from pyxsim.utils import parse_value, mylog
from astropy.table import Table
from yt.convenience import load
from yt.data_objects.api import Dataset

axis_wcs = [[1,2],[0,2],[0,1]]


def make_grid_source(ds, axis, width, center, redshift, area,
                     exp_time, source_model, sky_center, fov,
                     simput_prefix, depth=None, cosmology=None, dist=None,
                     absorb_model=None, nH=None, no_shifting=False,
                     sigma_pos=None, kernel="top_hat", overwrite=False,
                     prng=None):
    """
    Make a composite source from a rectangular grid of multiple
    fields of view. The photons from each field of view will be
    stored in a SIMPUT photon list, and all lists will be attached 
    to the same SIMPUT catalog.

    Parameters
    ----------
    ds : string or :class:`~yt.static_output.Dataset`
        The name of the dataset to generate the photons from, or
        an already loaded yt dataset.
    axis : string or integer
        The axis along which to project. "x", "y", or "z", or
        0, 1, or 2.
    width : tuple or a float.
        Width can have four different formats to support variable
        x and y widths.  They are:

        ==================================     =======================
        format                                 example
        ==================================     =======================
        (float, string)                        (10,'kpc')
        ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
        float                                  0.2
        (float, float)                         (0.2, 0.3)
        ==================================     =======================

        For example, (10, 'kpc') specifies a width that is 10 kiloparsecs
        wide in the x and y directions, ((10,'kpc'),(15,'kpc')) specifies a
        width that is 10 kiloparsecs wide along the x axis and 15
        kiloparsecs wide along the y axis.  In the other two examples, code
        units are assumed, for example (0.2, 0.3) specifies a width that has an
        x width of 0.2 and a y width of 0.3 in code units.
    center : string or array_like, optional
        The origin of the photon spatial coordinates. Accepts "c", "max", or
        a coordinate. If not specified, pyxsim attempts to use the "center"
        field parameter of the data_source.
    redshift : float
        The cosmological redshift for the photons.
    area : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The collecting area to determine the number of photons. If units are
        not specified, it is assumed to be in cm^2.
    exp_time : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The exposure time to determine the number of photons. If units are
        not specified, it is assumed to be in seconds.
    source_model : :class:`~pyxsim.source_models.SourceModel`
        A source model used to generate the photons.
    sky_center : array-like
        Center RA, Dec of the events in degrees.
    fov : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        Width of the field of view for each individual pointing of the
        grid. If units are specified, they are assumed to be in
        arcminutes.
    simput_prefix : string, optional
        The prefix of the SIMPUT catalog file to write or append 
        to. If not set, it will be the same as *prefix*.
    depth : A tuple or a float, optional
        A tuple containing the depth to project through and the string
        key of the unit: (width, 'unit').  If set to a float, code units
        are assumed
    dist : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`, optional
        The angular diameter distance, used for nearby sources. This may be
        optionally supplied instead of it being determined from the
        *redshift* and given *cosmology*. If units are not specified, it is
        assumed to be in kpc. To use this, the redshift must be set to zero.
    cosmology : :class:`~yt.utilities.cosmology.Cosmology`, optional
        Cosmological information. If not supplied, we try to get the
        cosmology from the dataset. Otherwise, LCDM with the default yt 
        parameters is assumed.
    absorb_model : string or :class:`~pyxsim.spectral_models.AbsorptionModel`
        A model for foreground galactic absorption, to simulate the
        absorption of events before being detected. This cannot be applied
        here if you already did this step previously in the creation of the
        :class:`~pyxsim.photon_list.PhotonList` instance. Known options for
        strings are "wabs" and "tbabs".
    nH : float, optional
        The foreground column density in units of 10^22 cm^{-2}. Only used
        if absorption is applied.
    no_shifting : boolean, optional
        If set, the photon energies will not be Doppler shifted.
    sigma_pos : float, optional
        Apply a gaussian smoothing operation to the sky positions of the
        events. This may be useful when the binned events appear blocky due
        to their uniform distribution within simulation cells. However, this
        will move the events away from their originating position on the
        sky, and so may distort surface brightness profiles and/or spectra.
        Should probably only be used for visualization purposes. Supply a
        float here to smooth with a standard deviation with this fraction
        of the cell size. Default: None
    kernel : string, optional
        The kernel used when smoothing positions of X-rays originating from
        SPH particles, "gaussian" or "top_hat". Default: "top_hat".
    overwrite : boolean, optional
        Set to True to overwrite a previous file.
    prng : integer or :class:`~numpy.random.RandomState` object 
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is to use the :mod:`numpy.random`
        module.
    """
    sky_center = np.array(sky_center)

    if not isinstance(ds, Dataset):
        ds = load(ds)

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

    mylog.info("Number of FOV to cover source is {:d} x {:d}.".format(nx, ny))

    axisx, axisy = axis_wcs[axis]

    first = True

    phlists = []
    ra = []
    dec = []
    for i in range(nx):
        for j in range(ny):
            box_center = center.copy()
            box_center[axisx] += (i-0.5*(nx-1))*fov_width
            box_center[axisy] += (j-0.5*(ny-1))*fov_width
            # Expand the box by a small amount to make sure we fill
            # the instrument FOV
            le = box_center - 0.5*fov_width*1.1
            re = box_center + 0.5*fov_width*1.1
            le[axis] = box_center[axis] - 0.5*depth
            re[axis] = box_center[axis] + 0.5*depth
            box = ds.box(le, re)
            photons = PhotonList.from_data_source(box, redshift, area,
                                                  exp_time, source_model,
                                                  center=box_center,
                                                  dist=dist, cosmology=cosmo)
            xsky = np.array([(box_center-center)[axisx]/D_A])
            ysky = np.array([(box_center-center)[axisy]/D_A])
            pixel_to_cel(xsky, ysky, sky_center)
            events = photons.project_photons("xyz"[axis], (xsky[0], ysky[0]),
                                             absorb_model=absorb_model,
                                             nH=nH, no_shifting=no_shifting,
                                             sigma_pos=sigma_pos, kernel=kernel,
                                             prng=prng)
            del photons

            phlist_prefix = "{:s}_{:d}_{:d}".format(simput_prefix, i, j)
            if first:
                append = False
                first = False
            else:
                append = True
            events.write_simput_file(phlist_prefix, overwrite=overwrite,
                                     append=append, simput_prefix=simput_prefix)

            del events

            phlists.append("{}_phlist.fits".format(phlist_prefix))
            ra.append(xsky[0])
            dec.append(ysky[0])

    t = Table([np.arange(nx*ny)+1, phlists, ra, dec],
              names=("src_id", "phlist", "ra", "dec"))
    t.meta["comments"] = ["simput: {}_simput.fits".format(simput_prefix),
                          "fov: {:.2f} arcmin".format(fov.v),
                          "sky_center: {:.2f},{:.2f}".format(sky_center[0], sky_center[1]),
                          "dims: {:d},{:d}".format(nx, ny)]
    t['ra'].format = '.2f'
    t['dec'].format = '.2f'

    outfile = "{}_photon_grid.txt".format(simput_prefix)
    mylog.info("Writing grid information to {}.".format(outfile))
    t.write(outfile, overwrite=overwrite, 
            delimiter="\t", format='ascii.commented_header')

    return outfile
