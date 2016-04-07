from photon_simulator import \
    PowerLawSourceModel, PhotonList, \
    XSpecAbsorbModel
from photon_simulator.tests.beta_model_source import \
    BetaModelSource
from yt.units.yt_array import YTQuantity
import numpy as np
from yt.testing import requires_file, requires_module
import os
import shutil
import tempfile
from yt.config import ytcfg
from yt.utilities.physical_constants import mp

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

xray_data_dir = ytcfg.get("yt", "xray_data_dir")

arf = os.path.join(xray_data_dir,"sxt-s_120210_ts02um_intallpxl.arf")
rmf = os.path.join(xray_data_dir,"ah_sxs_5ev_basefilt_20100712.rmf")

@requires_module("xspec")
@requires_file(arf)
@requires_file(rmf)
def test_power_law():
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

    bms = BetaModelSource()
    ds = bms.ds

    def _hard_emission(field, data):
        return YTQuantity(1.0e-27, "s**-1*keV**-1")*data["density"]/mp
    ds.add_field(("gas","hard_emission"), function=_hard_emission, units="keV**-1*s**-1*cm**-3")

    nH_sim = 0.1
    index_sim = 1.1
    abs_model = XSpecAbsorbModel("TBabs", nH_sim)

    A = YTQuantity(2000., "cm**2")
    exp_time = YTQuantity(2.0e5, "s")
    redshift = 0.01

    sphere = ds.sphere("c", (100.,"kpc"))

    plaw_model = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", index_sim, prng=bms.prng)

    photons = PhotonList.from_data_source(sphere, redshift, A, exp_time,
                                          plaw_model)

    D_A = photons.parameters["FiducialAngularDiameterDistance"]
    dist_fac = 1.0/(4.*np.pi*D_A*D_A*(1.+redshift)**2).in_cgs()
    norm_sim = float((sphere["hard_emission"]*sphere["cell_volume"].in_cgs()).sum()*dist_fac.in_cgs())

    events = photons.project_photons("z", responses=[arf,rmf],
                                     absorb_model=abs_model,
                                     convolve_energies=True, prng=bms.prng)
    events.write_spectrum("plaw_model_evt.pi", clobber=True)

    s = xspec.Spectrum("plaw_model_evt.pi")
    s.ignore("**-0.5")
    s.ignore("9.0-**")

    m = xspec.Model("tbabs*zpowerlw")
    m.zpowerlw.PhoIndex = 2.0
    m.zpowerlw.norm = 1.0
    m.zpowerlw.Redshift = redshift
    m.TBabs.nH = 0.02

    m.zpowerlw.PhoIndex.frozen = False
    m.zpowerlw.norm.frozen = False
    m.zpowerlw.Redshift.frozen = True
    m.TBabs.nH.frozen = True

    xspec.Fit.renorm()
    xspec.Fit.nIterations = 100
    xspec.Fit.perform()

    index = m.zpowerlw.PhoIndex.values[0]
    norm = m.zpowerlw.norm.values[0]

    dnorm = m.zpowerlw.norm.sigma
    dindex = m.zpowerlw.PhoIndex.sigma

    print(norm, norm_sim, dnorm)
    print(index, index_sim, dindex)

    assert np.abs(norm-norm_sim) < 1.645*dnorm
    assert np.abs(index-index_sim) < 1.645*dindex

    xspec.AllModels.clear()
    xspec.AllData.clear()

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_power_law()