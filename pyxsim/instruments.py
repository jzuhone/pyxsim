import numpy as np
from pyxsim.event_list import ConvolvedEventList
from pyxsim.utils import pyxsim_path
from soxs.instrument import \
    AuxiliaryResponseFile, \
    RedistributionMatrixFile
import os

aciss_arf = os.path.join(pyxsim_path, "response_files", "aciss_aimpt_cy18.arf")
aciss_rmf = os.path.join(pyxsim_path, "response_files", "aciss_aimpt_cy18.rmf")
sxs_arf = os.path.join(pyxsim_path, "response_files", "hitomi_sxs.arf")
sxs_rmf = os.path.join(pyxsim_path, "response_files", "hitomi_sxs.rmf")

class InstrumentSimulator(object):
    def __init__(self, name, arf_file, rmf_file):
        """
        Construct an instrument simulator.

        Parameters
        ----------
        inst_name : string
            The string corresponding to the name of the SOXS 
            instrument specification. 
        """
        self.name = name
        self.arf_file = arf_file
        self.rmf_file = rmf_file

    _arf = None
    @property
    def arf(self):
        if self._arf is None:
            self._arf = AuxiliaryResponseFile(self.arf_file)
        return self._arf

    _rmf = None
    @property
    def rmf(self):
        if self._rmf is None:
            self._rmf = RedistributionMatrixFile(self.rmf_file)
        return self._rmf

    def __call__(self, events, prng=None):
        if prng is None:
            prng = np.random
        flux = np.sum(events["eobs"]).to("erg") / \
               events.parameters["exp_time"]/events.parameters["area"]
        exp_time = events.parameters["exp_time"]
        emin = events["eobs"].min().value
        emax = events["eobs"].max().value
        new_events = self.arf.detect_events(events.events, exp_time, flux,
                                            [emin, emax], prng=prng)
        new_events = self.rmf.scatter_energies(new_events, prng=prng)
        parameters = {}
        parameters.update(events.parameters)
        parameters["channel_type"] = self.rmf.header["CHANTYPE"]
        parameters["mission"] = self.rmf.header.get("MISSION", "")
        parameters["instrument"] = self.rmf.header["INSTRUME"]
        parameters["telescope"] = self.rmf.header["TELESCOP"]
        parameters["arf"] = self.arf.filename
        parameters["rmf"] = self.rmf.filename
        return ConvolvedEventList(new_events, parameters, wcs=events.wcs)

# Specific instrument approximations

ACIS_S = InstrumentSimulator("acis-s", aciss_arf, aciss_rmf)
ACIS_I = InstrumentSimulator("acis-i", "acisi_aimpt_cy18.arf",
                             "acisi_aimpt_cy18.rmf")
Hitomi_SXS = InstrumentSimulator("hitomi_sxs", sxs_arf, sxs_rmf)
Athena_WFI = InstrumentSimulator("athena_wfi",
                                 "athena_wfi_1469_onaxis_w_filter_v20150326.arf",
                                 "athena_wfi_rmf_v20150326.rmf")
Athena_XIFU = InstrumentSimulator("athena_xifu",
                                  "athena_xifu_1469_onaxis_pitch249um_v20160401.arf",
                                  "athena_xifu_rmf_v20160401.rmf")
Lynx_Imager = InstrumentSimulator("lynx_hdxi", "xrs_hdxi_3x10.arf", "xrs_hdxi.rmf")
Lynx_Calorimeter = InstrumentSimulator("lynx_mucal", "xrs_mucal_3x10.arf", "xrs_mucal.rmf")
