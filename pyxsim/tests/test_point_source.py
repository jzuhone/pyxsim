from pyxsim.spectral_models import TBabsModel
from pyxsim.source_generators import make_point_sources
from yt.testing import requires_module
import os
import tempfile
import shutil
import numpy as np
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    get_fit_results
from soxs.spectra import Spectrum
from soxs.events import write_spectrum
from soxs.instrument_registry import get_instrument_from_registry, \
    make_simple_instrument
from soxs.instrument import RedistributionMatrixFile, \
    AuxiliaryResponseFile, instrument_simulator

prng = 49

try:
    make_simple_instrument("aciss_cy19", "sq_aciss_cy19", 20.0, 2400)
except KeyError:
    pass

acis_spec = get_instrument_from_registry("sq_aciss_cy19")

rmf = RedistributionMatrixFile(acis_spec["rmf"])
arf = AuxiliaryResponseFile(acis_spec['arf'])

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

    spec.apply_foreground_absorption(nH_sim, model="tbabs")

    positions = [(30.01, 45.0)]

    events = make_point_sources(area, exp_time, positions, (30.0, 45.0),
                                spec, prng=prng)

    events.write_simput_file("ptsrc", overwrite=True)

    instrument_simulator("ptsrc_simput.fits", "ptsrc_evt.fits",
                         exp_time, "sq_aciss_cy19", [30.0, 45.0],
                         overwrite=True, foreground=False, ptsrc_bkgnd=False,
                         instr_bkgnd=False,
                         prng=prng)

    write_spectrum("ptsrc_evt.fits", "point_source_evt.pi", overwrite=True)

    os.system("cp %s %s ." % (arf.filename, rmf.filename))

    load_user_model(mymodel, "tplaw")
    add_user_pars("tplaw", ["nH", "norm", "redshift", "alpha"],
                  [0.02, norm_sim*0.8, redshift, 0.9],
                  parmins=[0.0, 0.0, 0.0, 0.1],
                  parmaxs=[10.0, 1.0e9, 10.0, 10.0],
                  parfrozen=[True, False, True, False])

    load_pha("point_source_evt.pi")
    set_stat("cstat")
    set_method("simplex")
    ignore(":0.4, 9.0:")
    set_model("tplaw")
    fit()
    res = get_fit_results()

    assert np.abs(res.parvals[0]-norm_sim)/norm_sim < 0.05
    assert np.abs(res.parvals[1]-alpha_sim)/alpha_sim < 0.05

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_point_source()
