from pyxsim import \
    PowerLawSourceModel, make_photons, \
    WabsModel, project_photons, EventList
from pyxsim.tests.utils import \
    BetaModelSource
from yt.units.yt_array import YTQuantity
import numpy as np
import os
import shutil
import tempfile
from yt.utilities.cosmology import Cosmology 
from yt.utilities.physical_constants import mp
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    get_fit_results 
from soxs.instrument import RedistributionMatrixFile, \
    AuxiliaryResponseFile, instrument_simulator
from soxs.events import write_spectrum
from soxs.instrument_registry import get_instrument_from_registry, \
    make_simple_instrument


try:
    make_simple_instrument("chandra_acisi_cy0", "sq_acisi_cy0", 20.0, 2400)
except KeyError:
    pass

acis_spec = get_instrument_from_registry("sq_acisi_cy0")

rmf = RedistributionMatrixFile(acis_spec["rmf"])
arf = AuxiliaryResponseFile(acis_spec['arf'])


def mymodel(pars, x, xhi=None):
    dx = x[1]-x[0]
    xmid = x+0.5*dx
    wm = WabsModel(pars[0])
    wabs = wm.get_absorb(xmid)
    plaw = pars[1]*dx*(xmid*(1.0+pars[2]))**(-pars[3])
    return wabs*plaw


def test_power_law():
    plaw_fit(1.1, prng=33)
    plaw_fit(0.8, prng=28)
    plaw_fit(1.0, prng=23)


def plaw_fit(alpha_sim, prng=None):

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    bms = BetaModelSource()
    ds = bms.ds

    if prng is None:
        prng = bms.prng

    def _hard_emission(field, data):
        return data.ds.quan(1.0e-18, "s**-1*keV**-1")*data["density"]*data["cell_volume"]/mp
    ds.add_field(("gas", "hard_emission"), function=_hard_emission, 
                 units="keV**-1*s**-1", sampling_type="local")

    nH_sim = 0.02

    A = YTQuantity(2000., "cm**2")
    exp_time = YTQuantity(2.0e5, "s")
    redshift = 0.01

    cosmo = Cosmology()

    sphere = ds.sphere("c", (100., "kpc"))

    plaw_model = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", 
                                     alpha_sim, prng=prng)

    n_photons, n_cells = make_photons("plaw_photons.h5", sphere, redshift, 
                                      A, exp_time, plaw_model)

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")
    dist_fac = 1.0/(4.*np.pi*D_A*D_A*(1.+redshift)**3)
    norm_sim = float((sphere["hard_emission"]).sum()*dist_fac)*(1.+redshift)

    n_events = project_photons("plaw_photons.h5", "plaw_events.h5", "z", 
                               [30., 45.], absorb_model="wabs",
                               nH=nH_sim, prng=bms.prng, no_shifting=True)

    events = EventList("plaw_events.h5")

    events.write_to_simput("plaw", overwrite=True)

    instrument_simulator("plaw_simput.fits", "plaw_evt.fits",
                         exp_time, "sq_acisi_cy0", [30.0, 45.0],
                         overwrite=True, foreground=False, ptsrc_bkgnd=False,
                         instr_bkgnd=False,
                         prng=prng)

    write_spectrum("plaw_evt.fits", "plaw_model_evt.pi", overwrite=True)

    os.system("cp %s %s ." % (arf.filename, rmf.filename))

    load_user_model(mymodel, "wplaw")
    add_user_pars("wplaw", ["nH", "norm", "redshift", "alpha"],
                  [0.01, norm_sim*1.1, redshift, 0.9], 
                  parmins=[0.0, 0.0, 0.0, 0.1],
                  parmaxs=[10.0, 1.0e9, 10.0, 10.0],
                  parfrozen=[False, False, True, False])

    load_pha("plaw_model_evt.pi")
    set_stat("cstat")
    set_method("simplex")
    ignore(":0.6, 7.0:")
    set_model("wplaw")
    fit()
    res = get_fit_results()

    assert np.abs(res.parvals[0]-nH_sim)/nH_sim < 0.2
    assert np.abs(res.parvals[1]-norm_sim)/norm_sim < 0.05
    assert np.abs(res.parvals[2]-alpha_sim)/alpha_sim < 0.05

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


if __name__ == "__main__":
    test_power_law()
