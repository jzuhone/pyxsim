from yt.utilities.answer_testing.framework import \
    requires_ds, data_dir_load
from pyxsim.source_generators.xray_binaries import \
    make_xrb_photons, make_xrb_particles, \
    bolometric_correction
from yt.units.yt_array import uconcatenate
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

galaxy = "sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art"

@requires_ds(galaxy)
def test_xray_binaries():

    ds = data_dir_load(galaxy)

    sp = ds.sphere("max", (0.25, "Mpc"))

    metallicity_field = ("STAR", "metallicity_snii")
    scale_length = (1.0, "kpc")
    age_field = ("STAR", "age")

    new_ds = make_xrb_particles(sp, metallicity_field, age_field, 
                                scale_length=scale_length)

    dd = new_ds.all_data()

    area = ds.quan(25000.0, "cm**2")
    exp_time = ds.quan(300.0, "ks")
    emin = 0.1
    emax = 10.0

    photons_xrb = make_xrb_photons(new_ds, ds.current_redshift, area, 
                                   exp_time, emin, emax,
                                   center=sp.center, cosmology=ds.cosmology)

    D_L = ds.cosmology.luminosity_distance(0.0, ds.current_redshift)

    E = uconcatenate(photons_xrb["energy"])
    fluxes = dd["particle_luminosity"]*bolometric_correction/(4.0*np.pi*D_L*D_L)

    idxs = E*(1.0+ds.current_redshift) > 2.0
    E_mean = E[idxs].mean().to("erg")
    n1 = fluxes.sum()*exp_time*area/E_mean
    n2 = idxs.sum()
    dn = np.sqrt(n2)
    print(n1, n2, dn)
    assert np.abs(n1-n2) < 1.645*dn