from pyxsim.event_list import EventList
from soxs.spatial import FillFOVModel
from soxs.utils import parse_prng
from pyxsim.utils import parse_value
from yt.units.yt_array import YTArray
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    communication_system

comm = communication_system.communicators[-1]

def make_background(area, exp_time, fov, sky_center, spectrum, prng=None):
    r"""
    Create a new :class:`~pyxsim.event_list.EventList` which is filled
    uniformly with background events. 

    Parameters
    ----------
    area : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The collecting area to determine the number of events. If units are
        not specified, it is assumed to be in cm^2.
    exp_time : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The exposure time to determine the number of events. If units are
        not specified, it is assumed to be in seconds.
    fov : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The field of view of the event file. If units are not 
        provided, they are assumed to be in arcminutes.
    sky_center : array-like
        Center RA, Dec of the events in degrees.
    spectrum : :class:`~soxs.spectra.Spectrum`
        The spectrum for the background.
    prng : integer or :class:`~numpy.random.RandomState` object 
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is to use the :mod:`numpy.random` module.
    """
    prng = parse_prng(prng)

    fov = parse_value(fov, "arcmin")
    exp_time = parse_value(exp_time, "s")
    area = parse_value(area, "cm**2")

    t_exp = exp_time.value/comm.size

    e = spectrum.generate_energies(t_exp, area.value, prng=prng)
    fov_model = FillFOVModel(sky_center[0], sky_center[1], fov.value)
    ra, dec = fov_model.generate_coords(e.size, prng=prng)

    parameters = {"sky_center": YTArray(sky_center, "degree"),
                  "exp_time": exp_time,
                  "area": area}

    events = {}
    events["xsky"] = YTArray(ra.value, "degree")
    events["ysky"] = YTArray(dec.value, "degree")
    events["eobs"] = YTArray(e.value, "keV")

    return EventList(events, parameters)
