from pyxsim.event_list import EventList
from pyxsim.tests.utils import create_dummy_wcs
from pyxsim.spectral_models import XSpecAbsorbModel, XSpecThermalModel
from pyxsim.responses import AuxiliaryResponseFile, RedistributionMatrixFile
from yt.testing import requires_file, requires_module
import os
from yt.config import ytcfg
from numpy.random import RandomState
from yt.units.yt_array import YTQuantity
import tempfile
import shutil
import numpy as np

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
def test_background():

    import xspec

    xspec.Fit.statMethod = "cstat"
    xspec.Xset.addModelString("APECTHERMAL","yes")
    xspec.Fit.query = "yes"
    xspec.Fit.method = ["leven","10","0.01"]
    xspec.Fit.delta = 0.01
    xspec.Xset.chatter = 5

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    kT_sim = 1.0
    Z_sim = 0.0
    norm_sim = 4.0e-3

    exp_time = (200., "ks")

    arf = AuxiliaryResponseFile(arf_fn)
    rmf = RedistributionMatrixFile(rmf_fn)
    area = arf.max_area

    wcs = create_dummy_wcs()

    apec_model = XSpecThermalModel("apec", 0.1, 11.5, 20000,
                                   thermal_broad=False)
    abs_model = XSpecAbsorbModel("TBabs", 0.02)

    apec_model.prepare_spectrum(0.05)
    cspec, mspec = apec_model.get_spectrum(kT_sim)
    spec = (cspec+Z_sim*mspec)*YTQuantity(norm_sim*1.0e14, "cm**-5")
    ebins = apec_model.ebins

    params = {"ARF": arf_fn,
              "RMF": rmf_fn}
    events = EventList.create_empty_list(exp_time, area, wcs, parameters=params)

    new_events = events.add_background(ebins, spec, prng=prng,
                                       absorb_model=abs_model)

    new_events.convolve_energies(rmf, prng=prng)

    new_events.write_spectrum("background_evt.pi", clobber=True)

    s = xspec.Spectrum("background_evt.pi")
    s.ignore("**-0.5")
    s.ignore("9.0-**")

    m = xspec.Model("tbabs*apec")
    m.apec.kT = 2.0
    m.apec.Abundanc = 0.1
    m.apec.norm = 1.0
    m.apec.Redshift = 0.05
    m.TBabs.nH = 0.02

    m.apec.Abundanc.frozen = False
    m.apec.Redshift.frozen = True
    m.TBabs.nH.frozen = True

    xspec.Fit.renorm()
    xspec.Fit.nIterations = 100
    xspec.Fit.perform()

    kT  = m.apec.kT.values[0]
    Z = m.apec.Abundanc.values[0]
    norm = m.apec.norm.values[0]

    dkT = m.apec.kT.sigma
    dZ = m.apec.Abundanc.sigma
    dnorm = m.apec.norm.sigma

    assert np.abs(kT-kT_sim) < 1.645*dkT
    assert np.abs(Z-Z_sim) < 1.645*dZ
    assert np.abs(norm-norm_sim) < 1.645*dnorm

    print kT, kT_sim, dkT
    print Z, Z_sim, dZ
    print norm, norm_sim, dnorm

    xspec.AllModels.clear()
    xspec.AllData.clear()

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_background()