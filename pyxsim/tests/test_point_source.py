from pyxsim.event_list import EventList
from pyxsim.tests.utils import create_dummy_wcs
from pyxsim.instruments import ACIS_S
from pyxsim.spectral_models import TBabsModel
from yt.testing import requires_module
import os
from numpy.random import RandomState
import tempfile
import shutil
import numpy as np
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    covar, get_covar_results, set_covar_opt
from pyxsim.instruments import specs
from soxs.utils import write_spectrum
from numpy.testing import assert_allclose

prng = RandomState(24)

sigma_to_fwhm = 2.*np.sqrt(2.*np.log(2.))

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def mymodel(pars, x, xhi=None):
    dx = x[1]-x[0]
    xmid = x+0.5*dx
    wm = TBabsModel(pars[0])
    wabs = wm.get_absorb(xmid)
    plaw = pars[1]*dx*(xmid*(1.0+pars[2]))**(-pars[3])
    return wabs*plaw

@requires_module("sherpa")
def test_point_source():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    nH_sim = 0.02
    norm_sim = 1.0e-2
    alpha_sim = 0.95
    redshift = 0.02

    exp_time = (100., "ks")
    area = (3000., "cm**2")

    wcs = create_dummy_wcs()

    ebins = np.linspace(0.1, 11.5, 5001)
    emid = 0.5*(ebins[1:]+ebins[:-1])
    de = np.diff(ebins)
    spec = norm_sim*(emid*(1.0+redshift))**(-alpha_sim)*de

    abs_model = TBabsModel(nH_sim)

    events = EventList.create_empty_list(exp_time, area, wcs)

    positions = [(30.01, 45.0)]

    new_events = events.add_point_sources(positions, ebins, spec, prng=prng,
                                          absorb_model=abs_model)

    assert_allclose(np.unique(new_events["xsky"])[0].v, positions[0][0])
    assert_allclose(np.unique(new_events["ysky"])[0].v, positions[0][1])

    ACIS_S(new_events, "point_source_evt.fits", convolve_energies_only=True,
           instr_bkgnd=False, astro_bkgnd=False, prng=prng)

    os.system("cp %s ." % specs[ACIS_S.inst_name]["arf"])
    os.system("cp %s ." % specs[ACIS_S.inst_name]["rmf"])

    write_spectrum("point_source_evt.fits", "point_source_evt.pi", clobber=True)

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
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-nH_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-norm_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-alpha_sim) < res.parmaxes[2]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_point_source()
