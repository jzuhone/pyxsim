from pyxsim.event_list import EventList
from pyxsim.tests.utils import create_dummy_wcs
from pyxsim.spectral_models import TableApecModel, WabsModel
from pyxsim.instruments import ACIS_I
from yt.testing import requires_module
import os
from numpy.random import RandomState
import tempfile
import shutil
import numpy as np
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    covar, get_covar_results, set_covar_opt
from soxs.instrument import RedistributionMatrixFile
from soxs.events import write_spectrum
from pyxsim.instruments import specs

prng = RandomState(24)

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

rmf = RedistributionMatrixFile(specs[ACIS_I.inst_name]["rmf"])
fit_model = TableApecModel(rmf.elo[0], rmf.ehi[-1], rmf.n_de, thermal_broad=False)

def mymodel(pars, x, xhi=None):
    tm = WabsModel(pars[0])
    tbabs = tm.get_absorb(x)
    bapec = fit_model.return_spectrum(pars[1], pars[2], pars[3], pars[4])
    return tbabs*bapec

@requires_module("sherpa")
def test_background():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    kT_sim = 1.0
    Z_sim = 0.0
    norm_sim = 4.0e-2
    nH_sim = 0.04
    redshift = 0.01

    exp_time = (200., "ks")
    area = (1000., "cm**2")

    wcs = create_dummy_wcs()

    abs_model = WabsModel(nH_sim)

    events = EventList.create_empty_list(exp_time, area, wcs)

    spec_model = TableApecModel(0.05, 12.0, 5000, thermal_broad=False)
    spec = spec_model.return_spectrum(kT_sim, Z_sim, redshift, norm_sim)

    new_events = events.add_background(spec_model.ebins, spec, prng=prng,
                                       absorb_model=abs_model)

    ACIS_I(new_events, "background_evt.fits", convolve_energies_only=True,
           instr_bkgnd=False, foreground=False, ptsrc_bkgnd=False, prng=prng)

    os.system("cp %s ." % specs[ACIS_I.inst_name]["arf"])
    os.system("cp %s ." % specs[ACIS_I.inst_name]["rmf"])

    write_spectrum("background_evt.fits", "background_evt.pi", overwrite=True)

    load_user_model(mymodel, "wapec")
    add_user_pars("wapec", ["nH", "kT", "metallicity", "redshift", "norm"],
                  [0.01, 4.0, 0.2, redshift, norm_sim*0.8],
                  parmins=[0.0, 0.1, 0.0, -20.0, 0.0],
                  parmaxs=[10.0, 20.0, 10.0, 20.0, 1.0e9],
                  parfrozen=[False, False, False, True, False])

    load_pha("background_evt.pi")
    set_stat("cstat")
    set_method("simplex")
    ignore(":0.5, 8.0:")
    set_model("wapec")
    fit()
    set_covar_opt("sigma", 1.6)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-nH_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-kT_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-Z_sim) < res.parmaxes[2]
    assert np.abs(res.parvals[3]-norm_sim) < res.parmaxes[3]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_background()
