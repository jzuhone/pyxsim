from pyxsim import \
    PowerLawSourceModel, PhotonList, \
    XSpecAbsorbModel, AuxiliaryResponseFile, \
    RedistributionMatrixFile
from pyxsim.tests.beta_model_source import \
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

arf_fn = os.path.join(xray_data_dir,"sxt-s_120210_ts02um_intallpxl.arf")
rmf_fn = os.path.join(xray_data_dir,"ah_sxs_5ev_basefilt_20100712.rmf")

@requires_module("xspec")
@requires_file(arf_fn)
@requires_file(rmf_fn)
def test_power_law():
    plaw_fit(1.1)
    plaw_fit(0.8)
    plaw_fit(1.0)

def plaw_fit(alpha_sim):
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
        return YTQuantity(1.0e-19, "s**-1*keV**-1")*data["density"]*data["cell_volume"]/mp
    ds.add_field(("gas","hard_emission"), function=_hard_emission, units="keV**-1*s**-1")

    nH_sim = 0.02
    abs_model = XSpecAbsorbModel("TBabs", nH_sim)

    A = YTQuantity(2000., "cm**2")
    exp_time = YTQuantity(2.0e5, "s")
    redshift = 0.01

    sphere = ds.sphere("c", (100.,"kpc"))

    plaw_model = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", alpha_sim, prng=bms.prng)

    photons = PhotonList.from_data_source(sphere, redshift, A, exp_time,
                                          plaw_model)

    D_A = photons.parameters["FiducialAngularDiameterDistance"]
    dist_fac = 1.0/(4.*np.pi*D_A*D_A*(1.+redshift)**2).in_cgs()
    norm_sim = float((sphere["hard_emission"]).sum()*dist_fac.in_cgs())*(1.+redshift)

    arf = AuxiliaryResponseFile(arf_fn, rmffile=rmf_fn)
    rmf = RedistributionMatrixFile(rmf_fn)

    events = photons.project_photons("z", absorb_model=abs_model,
                                     prng=bms.prng, area_new=arf,
                                     no_shifting=True)
    events.convolve_energies(rmf, prng=bms.prng)
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

    alpha = m.zpowerlw.PhoIndex.values[0]
    norm = m.zpowerlw.norm.values[0]

    dnorm = m.zpowerlw.norm.sigma
    dalpha = m.zpowerlw.PhoIndex.sigma

    assert np.abs(alpha-alpha_sim) < 1.645*dalpha
    assert np.abs(norm-norm_sim) < 1.645*dnorm

    xspec.AllModels.clear()
    xspec.AllData.clear()

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_power_law()