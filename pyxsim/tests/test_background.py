from pyxsim.source_generators.background import make_background
from pyxsim.spectral_models import TableApecModel, WabsModel
from yt.testing import requires_module
from soxs.spectra import ApecGenerator
import os
import tempfile
import shutil
import numpy as np
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    get_fit_results
from soxs.instrument import RedistributionMatrixFile, \
    AuxiliaryResponseFile, instrument_simulator
from soxs.events import write_spectrum
from soxs.instrument_registry import get_instrument_from_registry, \
    make_simple_instrument

prng = 39

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

try:
    make_simple_instrument("acisi_cy19", "sq_acisi_cy19", 20.0, 2400)
except KeyError:
    pass

acis_spec = get_instrument_from_registry("sq_acisi_cy19")

rmf = RedistributionMatrixFile(acis_spec["rmf"])
arf = AuxiliaryResponseFile(acis_spec['arf'])

fit_model = TableApecModel(rmf.elo[0], rmf.ehi[-1], rmf.n_e, thermal_broad=False)

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
    fov = (10.0, "arcmin")

    prng = 24

    agen = ApecGenerator(0.05, 12.0, 5000, broadening=False)
    spec = agen.get_spectrum(kT_sim, Z_sim, redshift, norm_sim)
    spec.apply_foreground_absorption(norm_sim)

    events = make_background(area, exp_time, fov, (30.0, 45.0), spec, prng=prng)
    events.write_simput_file("bkgnd", overwrite=True)

    instrument_simulator("bkgnd_simput.fits", "bkgnd_evt.fits", 
                         exp_time, "sq_acisi_cy19", [30.0, 45.0],
                         overwrite=True, foreground=False, ptsrc_bkgnd=False,
                         instr_bkgnd=False,
                         prng=prng)

    write_spectrum("bkgnd_evt.fits", "background_evt.pi", overwrite=True)

    os.system("cp %s %s ." % (arf.filename, rmf.filename))

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
    res = get_fit_results()

    assert np.abs(res.parvals[0]-nH_sim)/nH_sim < 0.1
    assert np.abs(res.parvals[1]-kT_sim)/kT_sim < 0.05
    assert np.abs(res.parvals[2]-Z_sim) < 0.05
    assert np.abs(res.parvals[3]-norm_sim)/norm_sim < 0.05

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_background()
