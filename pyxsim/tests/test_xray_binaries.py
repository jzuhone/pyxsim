from yt.utilities.answer_testing.framework import \
    requires_ds, data_dir_load
from pyxsim.source_generators.xray_binaries import \
    make_xrb_photons, make_xrb_particles, \
    bc_lmxb, bc_hmxb, convert_bands, alpha_lmxb, \
    emin_lmxb, emax_lmxb, emin_hmxb, emax_hmxb
from yt.units.yt_array import uconcatenate, loadtxt
import numpy as np
import tempfile
import os
import shutil

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

galaxy = "FIRE_M12i_ref11/snapshot_600.hdf5"

@requires_ds(galaxy)
def test_xray_binaries():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = 25

    ds = data_dir_load(galaxy)

    def _age(field, data):
        z_s = 1.0 / data["PartType4", "StellarFormationTime"] - 1.0
        age = data.ds.cosmology.t_from_z(0.0) - data.ds.cosmology.t_from_z(z_s)
        age.convert_to_units("Gyr")
        return age

    ds.add_field(("PartType4", "particle_age"), function=_age, units="Gyr", 
                 particle_type=True)

    sp = ds.sphere("max", (0.25, "Mpc"))

    scale_length = (1.0, "kpc")
    age_field = ("PartType4", "particle_age")

    new_ds = make_xrb_particles(sp, age_field, scale_length, prng=prng)

    dd = new_ds.all_data()

    area = ds.quan(25000.0, "cm**2")
    exp_time = ds.quan(300.0, "ks")
    emin = 0.1
    emax = 10.0
    redshift = 0.01

    photons_xrb = make_xrb_photons(new_ds, redshift, area, 
                                   exp_time, emin, emax,
                                   center=sp.center,
                                   cosmology=ds.cosmology, prng=prng)

    D_L = ds.cosmology.luminosity_distance(0.0, redshift)

    kappa = convert_bands(alpha_lmxb, emin_hmxb, emax_hmxb,
                          emin_lmxb, emax_lmxb)

    E = uconcatenate(photons_xrb["energy"])

    lmxb_idxs = dd["particle_spectral_index"] < 2.0
    hmxb_idxs = np.logical_not(lmxb_idxs)

    flux = dd["particle_luminosity"][lmxb_idxs].sum()*kappa/bc_lmxb
    flux += dd["particle_luminosity"][hmxb_idxs].sum()/bc_hmxb
    flux /= 4.0*np.pi*D_L*D_L

    idxs = E*(1.0+redshift) > 2.0
    E_mean = E[idxs].mean().to("erg")
    n1 = flux*exp_time*area/E_mean
    n2 = idxs.sum()
    dn = np.sqrt(n2)

    assert np.abs(n1-n2) < 1.645*dn

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
