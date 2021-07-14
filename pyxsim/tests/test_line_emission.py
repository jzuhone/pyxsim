from pyxsim import \
    LineSourceModel, make_photons
from pyxsim.tests.utils import \
    BetaModelSource
from yt.units.yt_array import YTQuantity, YTArray
import numpy as np
from yt.utilities.cosmology import Cosmology
from yt.utilities.physical_constants import clight
import h5py
import os
import tempfile
import shutil

cross_section = YTQuantity(500.0e-22, "cm**3/s")
m_chi = YTQuantity(10.0, "GeV").to_equivalent("g", "mass_energy")


def test_line_emission():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    cosmo = Cosmology()

    bms = BetaModelSource()
    ds = bms.ds

    def _dm_emission(field, data):
        return (data["dark_matter_density"]/m_chi)**2*data["cell_volume"]*cross_section
    ds.add_field(("gas","dm_emission"), function=_dm_emission, units="s**-1",
                 sampling_type="cell")

    location = YTQuantity(3.5, "keV")
    sigma = YTQuantity(1000., "km/s")
    sigma_E = (location*sigma/clight).in_units("keV")

    A = YTQuantity(1000., "cm**2")
    exp_time = YTQuantity(2.0e5, "s")
    redshift = 0.01

    sphere = ds.sphere("c", (100.,"kpc"))

    line_model = LineSourceModel(location, "dm_emission", 
                                 sigma="dark_matter_dispersion", prng=32)

    n_photons, n_cells = make_photons("my_photons.h5", sphere, redshift, A, 
                                      exp_time, line_model)

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to("cm")

    dist_fac = 1.0/(4.*np.pi*D_A*D_A*(1.+redshift)**3)
    dm_E = (sphere["dm_emission"]).sum()

    with h5py.File("my_photons.h5", "r") as f:
        E = YTArray(f["data"]["energy"][:], "keV")
        n_E = len(E)

    n_E_pred = (exp_time*A*dm_E*dist_fac).in_units("dimensionless")

    loc = location/(1.+redshift)
    sig = sigma_E/(1.+redshift)

    assert np.abs(loc-E.mean()) < 1.645*sig/np.sqrt(n_E)
    assert np.abs(E.std()**2-sig*sig) < 1.645*np.sqrt(2*(n_E-1))*sig**2/n_E
    assert np.abs(n_E-n_E_pred) < 1.645*np.sqrt(n_E)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


if __name__ == "__main__":
    test_line_emission()