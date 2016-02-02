"""
A unit test for the photon_simulator analysis module.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.analysis_modules.photon_simulator.api import \
    XSpecThermalModel, XSpecAbsorbModel, \
    ThermalPhotonModel, PhotonList
from yt.analysis_modules.photon_simulator.tests.beta_model_source import \
    BetaModelSource
from yt.config import ytcfg
from yt.testing import requires_file, requires_module
import numpy as np
from yt.utilities.physical_constants import clight
import os
import tempfile
import shutil

ckms = clight.in_units("km/s").v

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

xray_data_dir = ytcfg.get("yt", "xray_data_dir")

arf = os.path.join(xray_data_dir,"sxt-s_120210_ts02um_intallpxl.arf")
rmf = os.path.join(xray_data_dir,"ah_sxs_5ev_basefilt_20100712.rmf")

@requires_module("xspec")
@requires_file(arf)
@requires_file(rmf)
def test_beta_model():
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

    X_H = 0.76

    A = 3000.
    exp_time = 1.0e5
    redshift = 0.05
    nH_sim = 0.02

    apec_model = XSpecThermalModel("bapec", 0.1, 11.5, 20000,
                                   thermal_broad=True)
    abs_model = XSpecAbsorbModel("TBabs", nH_sim)

    sphere = ds.sphere("c", (0.5, "Mpc"))
    sphere.set_field_parameter("X_H", X_H)

    mu_sim = -bms.v_shift / 1.0e5
    sigma_sim = bms.v_width / 1.0e5

    kT_sim = bms.kT
    Z_sim = bms.Z

    thermal_model = ThermalPhotonModel(apec_model, Zmet=Z_sim, X_H=X_H,
                                       prng=bms.prng)
    photons = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                      thermal_model)

    D_A = photons.parameters["FiducialAngularDiameterDistance"]

    norm_sim = sphere.quantities.total_quantity("emission_measure")
    norm_sim *= 1.0e-14/(4*np.pi*D_A*D_A*(1.+redshift)*(1.+redshift))
    norm_sim = float(norm_sim.in_cgs())

    events = photons.project_photons("z", responses=[arf,rmf],
                                     absorb_model=abs_model,
                                     convolve_energies=True, prng=bms.prng)
    events.write_spectrum("beta_model_evt.pi", clobber=True)

    s = xspec.Spectrum("beta_model_evt.pi")
    s.ignore("**-0.5")
    s.ignore("9.0-**")

    m = xspec.Model("tbabs*bapec")
    m.bapec.kT = 5.5
    m.bapec.Abundanc = 0.25
    m.bapec.norm = 1.0
    m.bapec.Redshift = 0.05
    m.bapec.Velocity = 300.0
    m.TBabs.nH = 0.02

    m.bapec.Velocity.frozen = False
    m.bapec.Abundanc.frozen = False
    m.bapec.Redshift.frozen = False
    m.TBabs.nH.frozen = True

    xspec.Fit.renorm()
    xspec.Fit.nIterations = 100
    xspec.Fit.perform()

    kT  = m.bapec.kT.values[0]
    mu = (m.bapec.Redshift.values[0]-redshift)*ckms
    Z = m.bapec.Abundanc.values[0]
    sigma = m.bapec.Velocity.values[0]
    norm = m.bapec.norm.values[0]

    dkT = m.bapec.kT.sigma
    dmu = m.bapec.Redshift.sigma*ckms
    dZ = m.bapec.Abundanc.sigma
    dsigma = m.bapec.Velocity.sigma
    dnorm = m.bapec.norm.sigma

    assert np.abs(mu-mu_sim) < dmu
    assert np.abs(kT-kT_sim) < dkT
    assert np.abs(Z-Z_sim) < dZ
    assert np.abs(sigma-sigma_sim) < dsigma
    assert np.abs(norm-norm_sim) < dnorm

    xspec.AllModels.clear()
    xspec.AllData.clear()

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
