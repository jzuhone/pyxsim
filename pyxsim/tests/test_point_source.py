from pyxsim.instruments import ACIS_S
from pyxsim.spectral_models import TBabsModel
from pyxsim.source_generators import make_point_sources
from yt.testing import requires_module
import os
import tempfile
import shutil
import numpy as np
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    covar, get_covar_results, set_covar_opt
from soxs.spectra import Spectrum

prng = 25

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def mymodel(pars, x, xhi=None):
    tm = TBabsModel(pars[0])
    tbabs = tm.get_absorb(x)
    dx = x[1]-x[0]
    plaw = pars[1]*dx*(x*(1.0+pars[2]))**(-pars[3])
    return tbabs*plaw

@requires_module("sherpa")
def test_point_source():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    nH_sim = 0.02
    norm_sim = 1.0e-4
    alpha_sim = 0.95
    redshift = 0.02

    exp_time = (100., "ks")
    area = (3000., "cm**2")

    spec = Spectrum.from_powerlaw(alpha_sim, redshift, norm_sim, 
                                  emin=0.1, emax=11.5, nbins=2000)

    de = 11.4/2000

    spec.apply_foreground_absorption(nH_sim, model="tbabs")

    positions = [(30.01, 45.0)]

    events = make_point_sources(area, exp_time, positions, (30.0, 45.0),
                                spec, prng=prng)

    new_events = ACIS_S(events, prng=prng)

    new_events.write_channel_spectrum("point_source_evt.pi", overwrite=True)

    os.system("cp %s %s ." % (ACIS_S.arf.filename, ACIS_S.rmf.filename))

    load_user_model(mymodel, "tplaw")
    add_user_pars("tplaw", ["nH", "norm", "redshift", "alpha"],
                  [0.01, norm_sim*0.8, redshift, 0.9],
                  parmins=[0.0, 0.0, 0.0, 0.1],
                  parmaxs=[10.0, 1.0e9, 10.0, 10.0],
                  parfrozen=[False, False, True, False])

    load_pha("point_source_evt.pi")
    set_stat("cstat")
    set_method("simplex")
    ignore(":0.5, 9.0:")
    set_model("tplaw")
    fit()
    set_covar_opt("sigma", 1.6)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-nH_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-norm_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-alpha_sim) < res.parmaxes[2]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_point_source()
