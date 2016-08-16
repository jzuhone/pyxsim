import time

from yt.analysis_modules.cosmological_observation.api import LightCone
from pyxsim.photon_list import PhotonList

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
        
        
    def generate_photons(self, area, exp_time, source_model, 
                         parameters=None, velocity_fields=None):
        pass