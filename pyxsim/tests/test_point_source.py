from pyxsim.event_list import EventList
from pyxsim.tests.utils import create_dummy_wcs
from pyxsim.spectral_models import XSpecAbsorbModel, XSpecThermalModel
from pyxsim.responses import AuxiliaryResponseFile
from yt.testing import requires_file, requires_module
import os
from yt.config import ytcfg
from numpy.random import RandomState
from yt.units.yt_array import YTQuantity

prng = RandomState(24)

xray_data_dir = ytcfg.get("yt", "xray_data_dir")

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

arf_fn = os.path.join(xray_data_dir, "acisi_aimpt_cy17.arf")
rmf_fn = os.path.join(xray_data_dir, "acisi_aimpt_cy17.rmf")

@requires_module("xspec")
@requires_file(arf_fn)
@requires_file(rmf_fn)
def test_point_source():

    exp_time = (100., "ks")

    arf = AuxiliaryResponseFile(arf_fn)
    area = arf.max_area

    wcs = create_dummy_wcs()

    apec_model = XSpecThermalModel("bapec", 0.1, 11.5, 20000,
                                   thermal_broad=True)
    abs_model = XSpecAbsorbModel("TBabs", 0.02)

    apec_model.prepare_spectrum(0.05)
    cspec, mspec = apec_model.get_spectrum(6.0)
    spec = (cspec+0.5*mspec)*YTQuantity(1.0e10, "cm**-5")
    ebins = apec_model.ebins

    params = {"ARF": arf_fn,
              "RMF": rmf_fn}
    events = EventList.create_empty_list(exp_time, area, wcs, parameters=params)

    positions = [(30.01, 45.0)]

    events.add_point_sources(positions, ebins, spec, prng=prng,
                             absorb_model=abs_model)

    events.write_fits_image("point_source.fits", clobber=True)

if __name__ == "__main__":
    test_point_source()