import time

from pyxsim.photon_list import PhotonList
from pyxsim.event_list import EventList
from pyxsim.utils import parse_value

try:
    from yt_astro_analysis.cosmological_observation.api import LightCone
except ImportError:
    try:
        from yt.analysis_modules.cosmological_observation.api import LightCone
    except ImportError:
        raise ImportError("Cannot import LightCone from yt or yt_astro_analysis!")

from yt.convenience import load
from yt.units.yt_array import uconcatenate, YTArray

from soxs.utils import parse_prng

from collections import defaultdict
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    communication_system

comm = communication_system.communicators[-1]

axes_lookup = [(1,2), (2,0), (0,1)]

class XrayLightCone(LightCone):
    def __init__(self, parameter_filename, simulation_type,
                 near_redshift, far_redshift, seed=None,
                 use_minimum_datasets=True, deltaz_min=0.0,
                 minimum_coherent_box_fraction=0.0):
        if seed is None:
            seed = time.time()
        super(XrayLightCone, self).__init__(parameter_filename, simulation_type,
                                            near_redshift, far_redshift,
                                            use_minimum_datasets=use_minimum_datasets,
                                            deltaz_min=deltaz_min,
                                            minimum_coherent_box_fraction=minimum_coherent_box_fraction)
        self.calculate_light_cone_solution(seed=seed)

    def generate_events(self, area, exp_time, angular_width,
                        source_model, sky_center, parameters=None,
                        velocity_fields=None, absorb_model=None,
                        nH=None, no_shifting=False, sigma_pos=None,
                        prng=None):
        """
        Generate projected events from a light cone simulation. 

        Parameters
        ----------
        area : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
            The collecting area to determine the number of events. If units are
            not specified, it is assumed to be in cm^2.
        exp_time : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
            The exposure time to determine the number of events. If units are
            not specified, it is assumed to be in seconds.
        angular_width : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
            The angular width of the light cone simulation. If units are not
            specified, it is assumed to be in degrees.
        source_model : :class:`~pyxsim.source_models.SourceModel`
            A source model used to generate the events.
        sky_center : array-like
            Center RA, Dec of the events in degrees.
        parameters : dict, optional
            A dictionary of parameters to be passed for the source model to use,
            if necessary.
        velocity_fields : list of fields
            The yt fields to use for the velocity. If not specified, the following will
            be assumed:
            ['velocity_x', 'velocity_y', 'velocity_z'] for grid datasets
            ['particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z'] for particle datasets
        absorb_model : string or :class:`~pyxsim.spectral_models.AbsorptionModel` 
            A model for foreground galactic absorption, to simulate the absorption
            of events before being detected. This cannot be applied here if you 
            already did this step previously in the creation of the 
            :class:`~pyxsim.photon_list.PhotonList` instance. Known options for 
            strings are "wabs" and "tbabs".
        nH : float, optional
            The foreground column density in units of 10^22 cm^{-2}. Only used if
            absorption is applied.
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
        prng : integer or :class:`~numpy.random.RandomState` object
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is to use the :mod:`numpy.random` module.
        """
        prng = parse_prng(prng)

        area = parse_value(area, "cm**2")
        exp_time = parse_value(exp_time, "s")
        aw = parse_value(angular_width, "deg")

        tot_events = defaultdict(list)

        for output in self.light_cone_solution:
            ds = load(output["filename"])
            ax = output["projection_axis"]
            c = output["projection_center"]*ds.domain_width + ds.domain_left_edge
            le = c.copy()
            re = c.copy()
            width = ds.quan(aw*output["box_width_per_angle"], "unitary").to("code_length")
            depth = ds.domain_width[ax].in_units("code_length")*output["box_depth_fraction"]
            le[ax] -= 0.5*depth
            re[ax] += 0.5*depth
            for off_ax in axes_lookup[ax]:
                le[off_ax] -= 0.5*width
                re[off_ax] += 0.5*width
            reg = ds.box(le, re)
            photons = PhotonList.from_data_source(reg, output['redshift'], area,
                                                  exp_time, source_model,
                                                  parameters=parameters,
                                                  center=c,
                                                  velocity_fields=velocity_fields,
                                                  cosmology=ds.cosmology)
            if sum(photons["num_photons"]) > 0:
                events = photons.project_photons("xyz"[ax], sky_center,
                                                 absorb_model=absorb_model, nH=nH,
                                                 no_shifting=no_shifting, 
                                                 sigma_pos=sigma_pos,
                                                 prng=prng)
                if events.num_events > 0:
                    tot_events["xsky"].append(events["xsky"])
                    tot_events["ysky"].append(events["ysky"])
                    tot_events["eobs"].append(events["eobs"])
                del events

            del photons

        parameters = {"exp_time": exp_time,
                      "area": area, 
                      "sky_center": YTArray(sky_center, "deg")}

        for key in tot_events:
            tot_events[key] = uconcatenate(tot_events[key])

        return EventList(tot_events, parameters)
