import time

from pyxsim.photon_list import PhotonList
from pyxsim.event_list import EventList
from pyxsim.utils import parse_value

from yt.analysis_modules.cosmological_observation.api import LightCone
from yt.convenience import load
from yt.units.yt_array import uconcatenate

from soxs.utils import parse_prng

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

        events_by_snapshot = []
        aw = parse_value(angular_width, "deg")

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
            if sum(photons["number_of_photons"]) > 0:
                events = photons.project_photons("xyz"[ax], sky_center,
                                                 absorb_model=absorb_model, nH=nH,
                                                 no_shifting=no_shifting, prng=prng)
                events_by_snapshot.append(events)
            del photons

        parameters = {}
        parameters.update(events_by_snapshot[-1].parameters)
        parameters.pop("d_a")
        parameters.pop("redshift")
        wcs = events_by_snapshot[-1].wcs

        tot_events = {}
        for key in ["eobs", "xsky", "ysky"]:
            tot_events[key] = events_by_snapshot[0][key]

        for events in events_by_snapshot[1:]:
            for key in ["eobs", "xsky", "ysky"]:
                tot_events[key] = uconcatenate(tot_events[key], events[key])

        x, y = wcs.wcs_world2pix(tot_events["xsky"].d, tot_events["ysky"].d, 1)
        tot_events["xpix"] = x
        tot_events["ypix"] = y

        return EventList(tot_events, parameters, wcs=wcs)
