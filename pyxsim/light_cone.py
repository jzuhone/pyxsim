import time

import astropy.wcs as pywcs

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
from yt.units.yt_array import uconcatenate, YTArray, YTQuantity

from soxs.utils import parse_prng

from collections import defaultdict

axes_lookup = [(1,2), (2,0), (0,1)]

class XrayLightCone(LightCone):
    def __init__(self, parameter_filename, simulation_type,
                 near_redshift, far_redshift, seed=None,
                 observer_redshift=0.0,
                 use_minimum_datasets=True, deltaz_min=0.0,
                 minimum_coherent_box_fraction=0.0):
        if seed is None:
            seed = time.time()
        super(XrayLightCone, self).__init__(parameter_filename, simulation_type,
                                            near_redshift, far_redshift,
                                            observer_redshift=observer_redshift,
                                            use_minimum_datasets=use_minimum_datasets,
                                            deltaz_min=deltaz_min,
                                            minimum_coherent_box_fraction=minimum_coherent_box_fraction)
        self.calculate_light_cone_solution(seed=seed)

    def generate_events(self, area, exp_time, angular_width, 
                        source_model, sky_center, parameters=None, 
                        velocity_fields=None, absorb_model=None, 
                        nH=None, no_shifting=False, prng=None):

        prng = parse_prng(prng)

        area = parse_value(area, "cm**2")
        exp_time = parse_value(exp_time, "s")
        aw = parse_value(angular_width, "deg")

        tot_events = defaultdict(list)

        dtheta = 1.0e88
        pix_center = None

        for output in self.light_cone_solution:
            ds = load(output["filename"])
            ax = output["projection_axis"]
            le = output["projection_center"].copy()
            re = output["projection_center"].copy()
            width = ds.quan(aw*output["box_width_per_angle"], "unitary").to("code_length")
            depth = ds.domain_width[ax].in_units("code_length")*output["box_depth_fraction"]
            le[ax] -= 0.5*depth
            re[ax] += 0.5*depth
            for off_ax in axes_lookup[ax]:
                le[off_ax] -= 0.5*width
                re[off_ax] += 0.5*width
            reg = ds.region(output["projection_center"], le, re)
            photons = PhotonList.from_data_source(reg, output['redshift'], area,
                                                  exp_time, source_model,
                                                  parameters=parameters,
                                                  velocity_fields=velocity_fields,
                                                  cosmology=ds.cosmology)
            if sum(photons["num_photons"]) > 0:
                events = photons.project_photons("xyz"[ax], sky_center,
                                                 absorb_model=absorb_model, nH=nH,
                                                 no_shifting=no_shifting, prng=prng)
                tot_events["xsky"].append(events["xsky"])
                tot_events["ysky"].append(events["ysky"])
                tot_events["eobs"].append(events["eobs"])
                if events.parameters["dtheta"].v < dtheta:
                    dtheta = events.parameters["dtheta"].v
                    pix_center = events.parameters["pix_center"]
                del events

            del photons

        wcs = pywcs.WCS(naxis=2)
        wcs.wcs.crpix = pix_center
        wcs.wcs.crval = list(sky_center)
        wcs.wcs.cdelt = [-dtheta, dtheta]
        wcs.wcs.ctype = ["RA---TAN","DEC--TAN"]
        wcs.wcs.cunit = ["deg"]*2

        parameters = {"exp_time": exp_time.v,
                      "area": area.v, 
                      "sky_center": YTArray(sky_center, "deg"),
                      "pix_center": pix_center,
                      "dtheta": YTQuantity(dtheta, "deg")}

        for key in tot_events:
            tot_events[key] = uconcatenate(tot_events[key])

        x, y = wcs.wcs_world2pix(tot_events["xsky"].d, tot_events["ysky"].d, 1)
        tot_events["xpix"] = x
        tot_events["ypix"] = y

        return EventList(tot_events, parameters, wcs=wcs)
