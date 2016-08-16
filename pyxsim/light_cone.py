import time

from pyxsim.photon_list import PhotonList
from pyxsim.utils import parse_value

from yt.analysis_modules.cosmological_observation.api import LightCone
from yt.convenience import load

axes_lookup = [(1,2),(2,0),(0,1)]

class XrayLightCone(LightCone):
    def __init__(self, parameter_filename, simulation_type,
                 near_redshift, far_redshift, seed=None,
                 observer_redshift=0.0,
                 use_minimum_datasets=True, deltaz_min=0.0,
                 minimum_coherent_box_fraction=0.0,
                 time_data=True, redshift_data=True,
                 find_outputs=False, set_parameters=None,
                 output_dir="LC", output_prefix="LightCone"):
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
                        source_model, parameters=None, 
                        velocity_fields=None, absorb_model=None, 
                        sky_center=None, no_shifting=False, 
                        prng=None):

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
            for a in axes_lookup[ax]:
                le[a] -= 0.5*width
                re[a] += 0.5*width
            reg = ds.region(output["projection_center"], le, re)
            photons = PhotonList.from_data_source(reg, output['redshift'], area,
                                                  exp_time, source_model,
                                                  parameters=parameters,
                                                  velocity_fields=velocity_fields,
                                                  cosmology=ds.cosmology)
            events = photons.project_photons("xyz"[ax], absorb_model=absorb_model, 
                                             sky_center=sky_center, 
                                             no_shifting=no_shifting, prng=prng)
            events_by_snapshot.append(events)

        tot_events = events_by_snapshot[0]
        for events in events_by_snapshot[1:]:
            tot_events = tot_events + events

        return tot_events
