from soxs.utils import parse_prng
from yt.funcs import ensure_list
from pyxsim.utils import parse_value
from yt.units.yt_array import uconcatenate, YTArray
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    communication_system
from pyxsim.event_list import EventList

comm = communication_system.communicators[-1]

def make_point_sources(area, exp_time, positions, sky_center,
                       spectra, prng=None):
    r"""
    Create a new :class:`~pyxsim.event_list.EventList` which contains
    point sources.

    Parameters
    ----------
    area : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The collecting area to determine the number of events. If units are
        not specified, it is assumed to be in cm^2.
    exp_time : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The exposure time to determine the number of events. If units are
        not specified, it is assumed to be in seconds.
    positions : array of source positions, shape 2xN
        The positions of the point sources in RA, Dec, where N is the
        number of point sources. Coordinates should be in degrees.
    sky_center : array-like
        Center RA, Dec of the events in degrees.
    spectra : list (size N) of :class:`~soxs.spectra.Spectrum` objects
        The spectra for the point sources, where N is the number 
        of point sources. Assumed to be in the observer frame.
    prng : integer or :class:`~numpy.random.RandomState` object 
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is to use the :mod:`numpy.random` module.
    """
    prng = parse_prng(prng)

    spectra = ensure_list(spectra)
    positions = ensure_list(positions)

    area = parse_value(area, "cm**2")
    exp_time = parse_value(exp_time, "s")

    t_exp = exp_time.value/comm.size

    x = []
    y = []
    e = []

    for pos, spectrum in zip(positions, spectra):
        eobs = spectrum.generate_energies(t_exp, area.value, prng=prng)
        ne = eobs.size
        x.append(YTArray([pos[0]] * ne, "deg"))
        y.append(YTArray([pos[1]] * ne, "deg"))
        e.append(YTArray.from_astropy(eobs))

    parameters = {"sky_center": YTArray(sky_center, "degree"),
                  "exp_time": exp_time,
                  "area": area}

    events = {}
    events["xsky"] = uconcatenate(x)
    events["ysky"] = uconcatenate(y)
    events["eobs"] = uconcatenate(e)

    return EventList(events, parameters)
