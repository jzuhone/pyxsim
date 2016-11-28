import numpy as np
from pyxsim.utils import pyxsim_path
from soxs.instrument import instrument_simulator
from soxs.instrument_registry import get_instrument_from_registry, \
    add_instrument_to_registry, instrument_registry
import uuid
import os

specs = {}
specs["acis-i"] = {"name": "acis-i",
                   "bkgnd": "acisi",
                   "fov": 67.1744,
                   "num_pixels": 8192,
                   "rmf": os.path.join(pyxsim_path, "response_files", "acisi_aimpt_cy18.rmf"),
                   "arf": os.path.join(pyxsim_path, "response_files", "acisi_aimpt_cy18.arf"),
                   "focal_length": 10.0,
                   "dither": True,
                   "psf": ["gaussian", 0.5]}

specs["acis-s"] = {"name": "acis-s",
                   "bkgnd": None,
                   "fov": 67.1744,
                   "num_pixels": 8192,
                   "rmf": os.path.join(pyxsim_path, "response_files", "aciss_aimpt_cy18.rmf"),
                   "arf": os.path.join(pyxsim_path, "response_files", "aciss_aimpt_cy18.arf"),
                   "focal_length": 10.0,
                   "dither": True,
                   "psf": ["gaussian", 0.5]}

specs["hitomi_sxs"] = {"name": "hitomi_sxs",
                       "psf": ["gaussian", 72.0],
                       "fov": 3.06450576,
                       "bkgnd": None,
                       "rmf": "",
                       "arf": "",
                       "focal_length": 5.6,
                       "dither": False,
                       "num_pixels": 6}

class InstrumentSimulator(object):
    def __init__(self, inst_name):
        """
        Construct an instrument simulator.

        Parameters
        ----------
        inst_name : string
            The string corresponding to the name of the SOXS 
            instrument specification. 
        """
        self.inst_name = inst_name

    def __call__(self, events, out_file, rebin=True,
                 convolve_psf=True, instr_bkgnd=True, 
                 astro_bkgnd=True, clobber=False, prng=None):
        if self.inst_name not in instrument_registry:
            add_instrument_to_registry(specs[self.inst_name])
        flux = np.sum(events["eobs"]).to("erg") / \
               events.parameters["ExposureTime"]/events.parameters["Area"]
        input_events = {"ra": events["xsky"].d, 
                        "dec": events["ysky"].d,
                        "energy": events["eobs"].d,
                        "flux": flux.v}
        inst = get_instrument_from_registry(self.inst_name)
        if prng is None:
            prng = np.random
        if not rebin:
            inst["num_pixels"] = int(2*events.parameters["pix_center"][0]-1.)
            inst["fov"] = events.parameters["dtheta"].v*60.0*inst["num_pixels"]
        if not convolve_psf:
            inst["psf"] = None
        inst["name"] = "_".join([inst["name"], uuid.uuid4().hex])
        add_instrument_to_registry(inst)
        exp_time = float(events.parameters["ExposureTime"])
        instrument_simulator(input_events, out_file, exp_time, inst["name"],
                             events.parameters["sky_center"].v, clobber=clobber, 
                             instr_bkgnd=instr_bkgnd, astro_bkgnd=astro_bkgnd, prng=prng)

# Specific instrument approximations

ACIS_S = InstrumentSimulator("acis-s")
ACIS_I = InstrumentSimulator("acis-i")
Hitomi_SXS = InstrumentSimulator("hitomi_sxs")
Athena_WFI = InstrumentSimulator("athena_wfi")
Athena_XIFU = InstrumentSimulator("athena_xifu")
XRS_Imager = InstrumentSimulator("hdxi")
XRS_Calorimeter = InstrumentSimulator("mucal")
