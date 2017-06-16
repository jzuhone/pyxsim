from pyxsim import \
    LineSourceModel, PhotonList
from pyxsim.tests.utils import \
    BetaModelSource
from yt.units.yt_array import YTQuantity, uconcatenate
import numpy as np
import yt.units as u
from yt.utilities.physical_constants import clight
from numpy.random import RandomState

cross_section = 500.0e-22*u.cm**3/u.s
m_chi = (10.0*u.GeV).to_equivalent("g", "mass_energy")

def test_line_emission():

    bms = BetaModelSource()
    ds = bms.ds

    def _dm_emission(field, data):
        return cross_section*(data["dark_matter_density"]/m_chi)**2*data["cell_volume"]
    ds.add_field(("gas","dm_emission"), function=_dm_emission, units="s**-1")

    location = YTQuantity(3.5, "keV")
    sigma = YTQuantity(1000., "km/s")
    sigma_E = (location*sigma/clight).in_units("keV")

    A = YTQuantity(1000., "cm**2")
    exp_time = YTQuantity(2.0e5, "s")
    redshift = 0.01

    sphere = ds.sphere("c", (100.,"kpc"))

    line_model = LineSourceModel(location, "dm_emission", 
                                 sigma="dark_matter_dispersion", prng=32)

    photons = PhotonList.from_data_source(sphere, redshift, A, exp_time,
                                          line_model)

    D_A = photons.parameters["fid_d_a"]
    dist_fac = 1.0/(4.*np.pi*D_A*D_A*(1.+redshift)**3)
    dm_E = (sphere["dm_emission"]).sum()

    E = uconcatenate(photons["energy"])
    n_E = len(E)

    n_E_pred = (exp_time*A*dm_E*dist_fac).in_units("dimensionless")

    loc = location/(1.+redshift)
    sig = sigma_E/(1.+redshift)

    assert np.abs(loc-E.mean()) < 1.645*sig/np.sqrt(n_E)
    assert np.abs(E.std()**2-sig*sig) < 1.645*np.sqrt(2*(n_E-1))*sig**2/n_E
    assert np.abs(n_E-n_E_pred) < 1.645*np.sqrt(n_E)

if __name__ == "__main__":
    test_line_emission()