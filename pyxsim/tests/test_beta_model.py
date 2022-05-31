"""
A unit test for the pyxsim analysis module.
"""

from pyxsim import \
    TableCIEModel, TBabsModel, \
    CIESourceModel, make_photons, \
    project_photons, EventList
from pyxsim.tests.utils import \
    BetaModelSource, ParticleBetaModelSource
import numpy as np
from yt.utilities.physical_constants import clight
import os
import tempfile
import shutil
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    get_fit_results
from soxs.instrument import RedistributionMatrixFile, \
    AuxiliaryResponseFile, instrument_simulator
from soxs.events import write_spectrum
from soxs.instrument_registry import get_instrument_from_registry, \
    make_simple_instrument
from yt.utilities.cosmology import Cosmology

cosmo = Cosmology()

ckms = clight.in_units("km/s").v

try:
    mucal_spec = get_instrument_from_registry("lynx_lxm")
    make_simple_instrument("lynx_lxm", "lynx_lxm_big", 20.0, 1200)
except KeyError:
    pass

rmf = RedistributionMatrixFile(mucal_spec["rmf"])
arf = AuxiliaryResponseFile(mucal_spec['arf'])
fit_model = TableCIEModel("apec", rmf.elo[0], rmf.ehi[-1], rmf.n_e)
agen_var = TableCIEModel("apec", rmf.elo[0], rmf.ehi[-1], rmf.n_e,
                         var_elem=["O", "Ca"], thermal_broad=True)


def mymodel(pars, x, xhi=None):
    dx = x[1]-x[0]
    tm = TBabsModel(pars[0])
    tbabs = tm.get_absorb(x+0.5*dx)
    bapec = fit_model.return_spectrum(pars[1], pars[2], pars[3], pars[4], velocity=pars[5])
    eidxs = np.logical_and(rmf.elo >= x[0]-0.5*dx, rmf.elo <= x[-1]+0.5*dx)
    return tbabs*bapec[eidxs]


def mymodel_var(pars, x, xhi=None):
    dx = x[1]-x[0]
    tm = TBabsModel(pars[0])
    tbabs = tm.get_absorb(x+0.5*dx)
    bapec = agen_var.return_spectrum(pars[1], pars[2], pars[3], pars[4],
                                     elem_abund={"O": pars[5], "Ca": pars[6]})
    eidxs = np.logical_and(rmf.elo >= x[0]-0.5*dx, rmf.elo <= x[-1]+0.5*dx)
    return tbabs*bapec[eidxs]


def test_beta_model():
    bms = BetaModelSource()
    do_beta_model(bms)


def test_beta_model_nomove():
    bms = BetaModelSource()
    do_beta_model(bms, axis="x", prng=89)


def test_beta_model_offaxis():
    bms = BetaModelSource()
    do_beta_model(bms, axis=[1.0, -2.0, 5.0], prng=78)


def test_particle_beta_model():
    bms = ParticleBetaModelSource()
    do_beta_model(bms, prng=29)


def test_particle_beta_model_nomove():
    bms = ParticleBetaModelSource()
    do_beta_model(bms, axis="x", prng=72)


def test_particle_beta_model_offaxis():
    bms = ParticleBetaModelSource()
    do_beta_model(bms, prng=67, axis=[1.0, -2.0, 5.0])


def do_beta_model(source, axis="z", prng=None):

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    if prng is None:
        prng = source.prng

    ds = source.ds

    A = 30000.
    exp_time = 1.0e4
    redshift = 0.05
    nH_sim = 0.02

    sphere = ds.sphere("c", (0.5, "Mpc"))

    kT_sim = source.kT
    Z_sim = source.Z

    thermal_model = CIESourceModel("apec", 0.1, 11.5, 20000, Z_sim, prng=prng)
    n_photons, n_cells = make_photons("my_photons", sphere, redshift, 
                                      A, exp_time, thermal_model)

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")

    norm_sim = sphere.quantities.total_quantity(("gas", "emission_measure"))
    norm_sim *= 1.0e-14/(4*np.pi*D_A*D_A*(1.+redshift)*(1.+redshift))
    norm_sim = float(norm_sim.in_cgs())

    v1, v2 = sphere.quantities.weighted_variance(("gas", "velocity_z"),
                                                 ("gas", "emission_measure"))

    if isinstance(axis, str):
        if axis == "z":
            fac = 1.0
        else:
            fac = 0.0
    else:
        axis /= np.sqrt(np.dot(axis, axis))
        fac = np.dot(axis, [0.0, 0.0, 1.0])

    sigma_sim = fac*float(v1.in_units("km/s"))
    mu_sim = -fac*float(v2.in_units("km/s"))

    n_events = project_photons("my_photons", "my_events", axis, [30.0, 45.0],
                               absorb_model="tbabs", nH=nH_sim, prng=prng)

    events = EventList("my_events.h5")

    events.write_to_simput("beta_model", overwrite=True)

    instrument_simulator("beta_model_simput.fits", "beta_model_evt.fits",
                         exp_time, "lynx_lxm_big", [30.0, 45.0],
                         overwrite=True, foreground=False, ptsrc_bkgnd=False,
                         instr_bkgnd=False, 
                         prng=prng)

    write_spectrum("beta_model_evt.fits", "beta_model_evt.pi", overwrite=True)

    os.system("cp %s %s ." % (arf.filename, rmf.filename))

    load_user_model(mymodel, "tbapec")
    add_user_pars("tbapec", ["nH", "kT", "metallicity", "redshift", "norm", "velocity"],
                  [0.02, 4.0, 0.2, 0.04, norm_sim*0.8, 300.0],
                  parmins=[0.0, 0.1, 0.0, -200.0, 0.0, 0.0],
                  parmaxs=[10.0, 20.0, 10.0, 200.0, 1.0e9, 20000.0],
                  parfrozen=[True, False, False, False, False, False])

    load_pha("beta_model_evt.pi")
    set_stat("cstat")
    set_method("levmar")
    ignore(":0.6, 8.0:")
    set_model("tbapec")
    fit()
    res = get_fit_results()

    redshift_sim = (1.0+mu_sim/ckms)*(1.0+redshift) - 1.0

    assert np.abs(res.parvals[0]-kT_sim)/kT_sim < 0.05
    assert np.abs(res.parvals[1]-Z_sim)/Z_sim < 0.05
    assert np.abs(res.parvals[2]-redshift_sim)/redshift_sim < 0.05
    assert np.abs(res.parvals[3]-norm_sim)/norm_sim < 0.05
    assert np.abs(res.parvals[4]-sigma_sim) < 30.0

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_vapec_beta_model():

    bms = BetaModelSource()

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = 47

    ds = bms.ds

    A = 30000.
    exp_time = 1.0e4
    redshift = 0.05
    nH_sim = 0.02

    sphere = ds.sphere("c", (0.5, "Mpc"))

    kT_sim = bms.kT
    Z_sim = bms.Z
    O_sim = bms.O
    Ca_sim = bms.Ca

    var_elem = {"O": ("stream", "oxygen"),
                "Ca": ("stream", "calcium")}

    thermal_model = CIESourceModel("apec", 0.1, 11.5, 20000, ("gas","metallicity"),
                                   var_elem=var_elem, prng=prng)

    n_photons = make_photons("my_photons", sphere, redshift, A, exp_time,
                             thermal_model)

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to("cm")

    norm_sim = sphere.quantities.total_quantity("emission_measure")
    norm_sim *= 1.0e-14/(4*np.pi*D_A*D_A*(1.+redshift)*(1.+redshift))
    norm_sim = float(norm_sim.in_cgs())

    n_events = project_photons("my_photons", "my_events", "z", [30.0, 45.0], 
                               absorb_model="tbabs", nH=nH_sim, prng=prng, 
                               no_shifting=True)

    events = EventList("my_events.h5")

    events.write_to_simput("vbeta_model", overwrite=True)

    instrument_simulator("vbeta_model_simput.fits", "vbeta_model_evt.fits",
                         exp_time, "lynx_lxm_big", [30.0, 45.0],
                         overwrite=True, foreground=False, ptsrc_bkgnd=False,
                         instr_bkgnd=False,
                         prng=prng)

    write_spectrum("vbeta_model_evt.fits", "vbeta_model_evt.pha", overwrite=True)

    os.system("cp %s %s ." % (arf.filename, rmf.filename))

    load_user_model(mymodel_var, "tbapec")
    add_user_pars("tbapec", ["nH", "kT", "abund", "redshift", "norm", "O", "Ca"],
                  [nH_sim, 4.0, Z_sim, redshift, norm_sim*0.8, 0.3, 0.5],
                  parmins=[0.0, 0.1, 0.0, -20.0, 0.0, 0.0, 0.0],
                  parmaxs=[10.0, 20.0, 10.0, 20.0, 1.0e9, 10.0, 10.0],
                  parfrozen=[True, False, True, True, False, False, False])

    load_pha("vbeta_model_evt.pha")
    set_stat("cstat")
    set_method("levmar")
    ignore(":0.6, 8.0:")
    set_model("tbapec")
    fit()
    res = get_fit_results()

    assert np.abs(res.parvals[0]-kT_sim)/kT_sim < 0.05
    assert np.abs(res.parvals[1]-norm_sim)/norm_sim < 0.05
    assert np.abs(res.parvals[2]-O_sim)/O_sim < 0.05
    assert np.abs(res.parvals[3]-Ca_sim)/Ca_sim < 0.05

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


if __name__ == "__main__":
    test_beta_model_nomove()
    test_beta_model_offaxis()
    test_particle_beta_model_nomove()
    test_particle_beta_model_offaxis()
    test_beta_model()
    test_particle_beta_model()
    test_vapec_beta_model()
