import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from six import string_types
from soxs.utils import parse_prng
from yt.loaders import load_particles
from yt.units.yt_array import YTArray, YTQuantity, uconcatenate
from yt.utilities.physical_ratios import keV_per_erg

from pyxsim.photon_list import make_photons
from pyxsim.source_models import PowerLawSourceModel
from pyxsim.utils import mylog, parse_value

"""
Papers used in this code:

Gilfanov, M. 2004, MNRAS, 349, 146

Mineo, S., Gilfanov, M., & Sunyaev, R. 2012, MNRAS, 419, 2095
"""


# Function to calculate the scale factor for a power
# law with F = K*E**-alpha (K in units of ct/s/keV)
def get_scale_factor(ind, emin, emax):
    if ind == 2.0:
        k = np.log(emax / emin)
    else:
        k = (emax ** (2.0 - ind) - emin ** (2.0 - ind)) / (2.0 - ind)
    return keV_per_erg / k


# Function to convert between two different energy
# bands for a single power law
def convert_bands(ind, emin_a, emax_a, emin_b, emax_b):
    if ind == 2.0:
        k = np.log(emax_a / emin_a)
        k /= np.log(emax_b / emin_b)
    else:
        k = emax_a ** (2.0 - ind) - emin_a ** (2.0 - ind)
        k /= emax_b ** (2.0 - ind) - emin_b ** (2.0 - ind)
    return k


# Spectral indices for both types of XRBs

alpha_lmxb = 1.56
alpha_hmxb = 2.0

# Energy bands for luminosities in XRB
# distribution functions

emin_lmxb = 0.5
emax_lmxb = 8.0

emin_hmxb = 2.0
emax_hmxb = 10.0

# Bolometric corrections

bc_lmxb = convert_bands(alpha_lmxb, 0.03, 100.0, emin_lmxb, emax_lmxb)
bc_hmxb = convert_bands(alpha_hmxb, 0.03, 100.0, emin_hmxb, emax_hmxb)

# Range of luminosities common to both types of XRBs

Lmin = 1.0e-3
Lcut = 1000.0
nbins = 1000
Lbins = np.logspace(np.log10(Lmin), np.log10(Lcut), nbins + 1)
logLbins = np.log10(Lbins)
logLmid = 0.5 * (logLbins[1:] + logLbins[:-1])

# LMXB distribution function from Gilfanov 2004

alpha1 = 1.0
alpha2 = 1.86
alpha3 = 4.8

# The luminosities from Gilfanov 2004 are
# in the 0.5-8 keV band.

Lb1 = 0.19
Lb2 = 5.0

K1 = 440.4
K2 = K1 * (Lb1 / Lb2) ** alpha2
K3 = K2 * (Lb2 / Lcut) ** alpha3

C1 = Lb1 * K1
C2 = Lb2 * K2 / (1.0 - alpha2)
C3 = Lcut * K3 / (1.0 - alpha3)

D2 = C2 * (Lb1 / Lb2) ** (1.0 - alpha2)
D3 = C3 * (Lb2 / Lcut) ** (1.0 - alpha3)

I1 = C1 * np.log(Lb1 / Lmin)
I2 = C2 - D2 + I1
I3 = C3 - D3 + I2


def lmxb_cdf(L):
    if L < Lb1:
        N = C1 * np.log(L / Lmin)
    elif Lb1 <= L < Lb2:
        N = C2 * (L / Lb2) ** (1.0 - alpha2) - D2 + I1
    elif Lb2 <= L < Lcut:
        N = C3 * (L / Lcut) ** (1.0 - alpha3) - D3 + I2
    else:
        N = I3
    return N


# HMXB distribution function from Mineo et al. 2012

chi = 1.49
gamma1 = 1.58
gamma2 = 2.73
Lb = 110.0
A = Lb ** (gamma2 - gamma1)

E1 = chi / (1.0 - gamma1)
E2 = chi * A / (1.0 - gamma2)

F1 = E1 * Lmin ** (1.0 - gamma1)
F2 = E2 * Lb ** (1.0 - gamma2)

J1 = E1 * Lb ** (1.0 - gamma1) - F1
J2 = E2 * Lcut ** (1.0 - gamma2) - F2 + J1


def hmxb_cdf(L):
    if L < Lb:
        N = E1 * L ** (1.0 - gamma1) - F1
    elif Lb <= L < Lcut:
        N = E2 * L ** (1.0 - gamma2) - F2 + J1
    else:
        N = J2
    return N


def make_xrb_particles(
    data_source, age_field, scale_length, sfr_time_range=(1.0, "Gyr"), prng=None
):
    r"""
    This routine generates an in-memory dataset composed of X-ray binary particles
    from an input data source containing star particles.

    Parameters
    ----------
    data_source : :class:`~yt.data_objects.data_containers.YTSelectionContainer`
        The yt data source to obtain the data from, such as a sphere, box, disk,
        etc.
    age_field : string or (type, name) field tuple
        The stellar age field. Must be in some kind of time units.
    scale_length : string, (ftype, fname) tuple, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The radial length scale over which to scatter the XRB particles
        from their parent star particle. Can be the name of a smoothing
        length field for the stars, a (value, unit) tuple, or a YTQuantity.
    sfr_time_range : string, (ftype, fname) tuple, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`, optional
        The recent time range over which to calculate the star formation rate from
        the current time in the dataset. Default: 1.0 Gyr
    prng : integer or :class:`~numpy.random.RandomState` object
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is to use the :mod:`numpy.random` module.
    """
    prng = parse_prng(prng)

    ds = data_source.ds

    ptype = data_source._determine_fields(age_field)[0][0]

    t = data_source[age_field].to("Gyr")
    m = data_source[(ptype, "particle_mass")].to("Msun")

    sfr_time_range = parse_value(sfr_time_range, "Gyr")

    recent = t < sfr_time_range

    n_recent = recent.sum()

    if n_recent == 0:
        sfr = 0.0
    else:
        sfr = (m[recent].sum() / sfr_time_range).to("Msun/yr").v

    mylog.info(
        "%d star particles were formed in the last %s for a SFR of %4.1f Msun/yr.",
        n_recent,
        sfr_time_range,
        sfr,
    )

    mtot = m.sum()

    npart = m.size

    scale_field = None
    if isinstance(scale_length, tuple):
        if isinstance(scale_length[0], string_types):
            scale_field = scale_length
    elif isinstance(scale_length, string_types):
        scale_field = (ptype, scale_length)

    if scale_field is None:
        if isinstance(scale_length, tuple):
            scale = YTArray([scale_length[0]] * npart, scale_length[1])
        elif isinstance(scale_length, YTQuantity):
            scale = YTArray([scale_length] * npart)
        else:
            scale = YTArray([scale_length[0]] * npart, "kpc")
    else:
        scale = data_source[scale_length]

    scale = scale.to("kpc").d

    N_l = lmxb_cdf(Lcut) * mtot.v * 1.0e-11
    N_h = hmxb_cdf(Lcut) * sfr

    N_all = N_l + N_h

    if N_all == 0.0:
        raise RuntimeError("There are no X-ray binaries to generate!")

    # Compute conversion factors from luminosity to count rate

    lmxb_factor = get_scale_factor(alpha_lmxb, emin_lmxb, emax_lmxb)
    hmxb_factor = get_scale_factor(alpha_hmxb, emin_hmxb, emax_hmxb)

    xp = []
    yp = []
    zp = []
    vxp = []
    vyp = []
    vzp = []
    lp = []
    rp = []
    ap = []

    if N_l > 0.0:

        F_l = np.zeros(nbins + 1)
        for i in range(1, nbins + 1):
            F_l[i] = lmxb_cdf(Lbins[i])
        F_l /= F_l[-1]
        invcdf_l = InterpolatedUnivariateSpline(F_l, logLbins)

        n_l = prng.poisson(lam=N_l * m / mtot)

        mylog.info("Number of low-mass X-ray binaries: %s", n_l.sum())

        for i, n in enumerate(n_l):
            if n > 0:
                randvec = prng.uniform(size=n)
                l = YTArray(10 ** invcdf_l(randvec) * 1.0e38, "erg/s")
                r = YTArray(l.v * lmxb_factor, "photons/s/keV")
                # Now convert output luminosities to bolometric
                l *= bc_lmxb
                x = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                y = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                z = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                x += data_source[ptype, "particle_position_x"][i].to("kpc")
                y += data_source[ptype, "particle_position_y"][i].to("kpc")
                z += data_source[ptype, "particle_position_z"][i].to("kpc")
                vx = YTArray([data_source[ptype, "particle_velocity_x"][i]] * n).to(
                    "km/s"
                )
                vy = YTArray([data_source[ptype, "particle_velocity_y"][i]] * n).to(
                    "km/s"
                )
                vz = YTArray([data_source[ptype, "particle_velocity_z"][i]] * n).to(
                    "km/s"
                )
                xp.append(x)
                yp.append(y)
                zp.append(z)
                vxp.append(vx)
                vyp.append(vy)
                vzp.append(vz)
                lp.append(l)
                rp.append(r)
                ap.append(np.array([alpha_lmxb] * n))

    if N_h > 0.0:

        F_h = np.zeros(nbins + 1)
        for i in range(1, nbins + 1):
            F_h[i] = hmxb_cdf(Lbins[i])
        F_h /= F_h[-1]
        invcdf_h = InterpolatedUnivariateSpline(F_h, logLbins)

        n_h = prng.poisson(lam=N_h * m / mtot)

        mylog.info("Number of high-mass X-ray binaries: %s", n_h.sum())

        for i, n in enumerate(n_h):
            if n > 0:
                randvec = prng.uniform(size=n)
                l = YTArray(10 ** invcdf_h(randvec) * 1.0e38, "erg/s")
                r = YTArray(l.v * hmxb_factor, "photons/s/keV")
                # Now convert output luminosities to bolometric
                l *= bc_hmxb
                x = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                y = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                z = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                x += data_source[ptype, "particle_position_x"][i].to("kpc")
                y += data_source[ptype, "particle_position_y"][i].to("kpc")
                z += data_source[ptype, "particle_position_z"][i].to("kpc")
                vx = YTArray([data_source[ptype, "particle_velocity_x"][i]] * n).to(
                    "km/s"
                )
                vy = YTArray([data_source[ptype, "particle_velocity_y"][i]] * n).to(
                    "km/s"
                )
                vz = YTArray([data_source[ptype, "particle_velocity_z"][i]] * n).to(
                    "km/s"
                )
                xp.append(x)
                yp.append(y)
                zp.append(z)
                vxp.append(vx)
                vyp.append(vy)
                vzp.append(vz)
                lp.append(l)
                rp.append(r)
                ap.append(np.array([alpha_hmxb] * n))

    xp = uconcatenate(xp)
    yp = uconcatenate(yp)
    zp = uconcatenate(zp)
    vxp = uconcatenate(vxp)
    vyp = uconcatenate(vyp)
    vzp = uconcatenate(vzp)
    lp = uconcatenate(lp)
    rp = uconcatenate(rp)
    ap = uconcatenate(ap)

    data = {
        "particle_position_x": (xp.d, str(xp.units)),
        "particle_position_y": (yp.d, str(yp.units)),
        "particle_position_z": (zp.d, str(zp.units)),
        "particle_velocity_x": (vxp.d, str(vxp.units)),
        "particle_velocity_y": (vyp.d, str(vyp.units)),
        "particle_velocity_z": (vzp.d, str(vzp.units)),
        "particle_luminosity": (lp.d, str(lp.units)),
        "particle_count_rate": (rp.d, str(rp.units)),
        "particle_spectral_index": ap,
    }

    dle = ds.domain_left_edge.to("kpc").v
    dre = ds.domain_right_edge.to("kpc").v

    bbox = np.array([[dle[i], dre[i]] for i in range(3)])

    new_ds = load_particles(
        data,
        bbox=bbox,
        length_unit="kpc",
        time_unit="Myr",
        mass_unit="Msun",
        velocity_unit="km/s",
    )

    return new_ds


def make_xrb_photons(
    photon_prefix,
    ds,
    redshift,
    area,
    exp_time,
    emin,
    emax,
    center="c",
    cosmology=None,
    prng=None,
):
    r"""
    Take a dataset produced by
    :func:`~pyxsim.source_generators.xray_binaries.make_xrb_particles`
    and produce a :class:`~pyxsim.photon_list.PhotonList`.

    Parameters
    ----------
    photon_prefix : string
        The prefix of the filename(s) containing the photon list. If run in
        serial, the filename will be "{photon_prefix}.h5", if run in
        parallel, the filenames will be "{photon_prefix}.{mpi_rank}.h5".
    ds : :class:`~yt.data_objects.static_output.Dataset`
        The dataset of XRB particles to use to make the photons.
    redshift : float
        The cosmological redshift for the photons.
    area : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The collecting area to determine the number of photons. If units are
        not specified, it is assumed to be in cm^2.
    exp_time : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The exposure time to determine the number of photons. If units are
        not specified, it is assumed to be in seconds.
    emin : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The minimum energy of the photons to be generated, in the rest frame of
        the source. If units are not given, they are assumed to be in keV.
    emax : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The maximum energy of the photons to be generated, in the rest frame of
        the source. If units are not given, they are assumed to be in keV.
    center : string or array_like, optional
        The origin of the photon spatial coordinates. Accepts "c", "max", or a coordinate.
        If not specified, pyxsim attempts to use the "center" field parameter of the data_source.
    cosmology : :class:`~yt.utilities.cosmology.Cosmology`, optional
        Cosmological information. If not supplied, we try to get
        the cosmology from the dataset. Otherwise, LCDM with
        the default yt parameters is assumed.
    prng : integer or :class:`~numpy.random.RandomState` object
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is to use the :mod:`numpy.random` module.
    """
    dd = ds.all_data()
    e0 = (1.0, "keV")
    prng = parse_prng(prng)
    xrb_model = PowerLawSourceModel(
        e0, emin, emax, "particle_count_rate", "particle_spectral_index", prng=prng
    )
    n_photons, n_cells = make_photons(
        photon_prefix,
        dd,
        redshift,
        area,
        exp_time,
        xrb_model,
        center=center,
        point_sources=True,
        cosmology=cosmology,
    )
    return n_photons, n_cells
