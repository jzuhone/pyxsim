"""
A unit test for the pyxsim analysis module.
"""

from pyxsim import \
    TableApecModel, TBabsModel, \
    ThermalSourceModel, PhotonList, \
    Lynx_Calorimeter
from pyxsim.tests.utils import \
    BetaModelSource, ParticleBetaModelSource
from yt.testing import requires_module
import numpy as np
from yt.utilities.physical_constants import clight
import os
import tempfile
import shutil
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    covar, get_covar_results, set_covar_opt, thaw
from soxs.utils import convert_rmf
from soxs.instrument import RedistributionMatrixFile, AuxiliaryResponseFile
from soxs.instrument_registry import get_instrument_from_registry

ckms = clight.in_units("km/s").v

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

mucal_spec = get_instrument_from_registry("mucal")

rmf = RedistributionMatrixFile(mucal_spec["rmf"])
arf = AuxiliaryResponseFile(mucal_spec['arf'])
fit_model = TableApecModel(rmf.elo[0], rmf.ehi[-1], rmf.n_de)
agen_var = TableApecModel(rmf.elo[0], rmf.ehi[-1], rmf.n_de,
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

@requires_module("sherpa")
def test_beta_model():
    bms = BetaModelSource()
    do_beta_model(bms, "velocity_z", "emission_measure")

@requires_module("sherpa")
def test_particle_beta_model():
    bms = ParticleBetaModelSource()
    do_beta_model(bms, "particle_velocity_z", 
                  ("io", "emission_measure"), prng=29)

def do_beta_model(source, v_field, em_field, prng=None):

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

    thermal_model = ThermalSourceModel("apec", 0.1, 11.5, 20000, 
                                       Zmet=Z_sim, prng=prng)
    photons = PhotonList.from_data_source(sphere, redshift, A, exp_time,
                                          thermal_model)

    D_A = photons.parameters["fid_d_a"]

    norm_sim = sphere.quantities.total_quantity(em_field)
    norm_sim *= 1.0e-14/(4*np.pi*D_A*D_A*(1.+redshift)*(1.+redshift))
    norm_sim = float(norm_sim.in_cgs())

    v1, v2 = sphere.quantities.weighted_variance(v_field, em_field)
    sigma_sim = float(v1.in_units("km/s"))
    mu_sim = -float(v2.in_units("km/s"))

    events = photons.project_photons("z", [30.0, 45.0], absorb_model="tbabs",
                                     nH=nH_sim, prng=prng)

    new_events = Lynx_Calorimeter(events, prng=prng)

    os.system("cp %s %s ." % (arf.filename, rmf.filename))
    convert_rmf(rmf.filename)

    new_events.write_channel_spectrum("beta_model_evt.pi", overwrite=True)

    load_user_model(mymodel, "tbapec")
    add_user_pars("tbapec", ["nH", "kT", "metallicity", "redshift", "norm", "velocity"],
                  [0.01, 4.0, 0.2, 0.04, norm_sim*0.8, 300.0],
                  parmins=[0.0, 0.1, 0.0, -200.0, 0.0, 0.0],
                  parmaxs=[10.0, 20.0, 10.0, 200.0, 1.0e9, 20000.0],
                  parfrozen=[False, False, False, False, False, False])

    load_pha("beta_model_evt.pi")
    set_stat("cstat")
    set_method("levmar")
    ignore(":0.6, 8.0:")
    set_model("tbapec")
    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    redshift_sim = (1.0+mu_sim/ckms)*(1.0+redshift) - 1.0

    assert np.abs(res.parvals[0]-nH_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-kT_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-Z_sim) < res.parmaxes[2]
    assert np.abs(res.parvals[3]-redshift_sim) < res.parmaxes[3]
    assert np.abs(res.parvals[4]-norm_sim) < res.parmaxes[4]
    assert np.abs(res.parvals[5]-sigma_sim) < res.parmaxes[5]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_vapec_beta_model():

    bms = BetaModelSource()

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = 45

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

    thermal_model = ThermalSourceModel("apec", 0.1, 11.5, 20000,
                                       var_elem=var_elem,
                                       Zmet=("gas","metallicity"), 
                                       prng=prng)

    photons = PhotonList.from_data_source(sphere, redshift, A, exp_time,
                                          thermal_model)

    D_A = photons.parameters["fid_d_a"]

    norm_sim = sphere.quantities.total_quantity("emission_measure")
    norm_sim *= 1.0e-14/(4*np.pi*D_A*D_A*(1.+redshift)*(1.+redshift))
    norm_sim = float(norm_sim.in_cgs())

    events = photons.project_photons("z", [30.0, 45.0], absorb_model="tbabs",
                                     nH=nH_sim, prng=prng, no_shifting=True)

    new_events = Lynx_Calorimeter(events, prng=prng)

    os.system("cp %s %s ." % (arf.filename, rmf.filename))
    convert_rmf(rmf.filename)

    new_events.write_channel_spectrum("var_abund_beta_model_evt.pha", overwrite=True)

    load_user_model(mymodel_var, "tbapec")
    add_user_pars("tbapec", ["nH", "kT", "abund", "redshift", "norm", "O", "Ca"],
                  [nH_sim, 4.0, Z_sim, redshift, norm_sim*0.8, 0.3, 0.5],
                  parmins=[0.0, 0.1, 0.0, -20.0, 0.0, 0.0, 0.0],
                  parmaxs=[10.0, 20.0, 10.0, 20.0, 1.0e9, 10.0, 10.0],
                  parfrozen=[True, False, True, True, False, False, False])

    load_pha("var_abund_beta_model_evt.pha")
    set_stat("cstat")
    set_method("levmar")
    ignore(":0.6, 8.0:")
    set_model("tbapec")
    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-kT_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-norm_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-O_sim) < res.parmaxes[2]
    assert np.abs(res.parvals[3]-Ca_sim) < res.parmaxes[3]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    #test_beta_model()
    #test_particle_beta_model()
    test_vapec_beta_model()
