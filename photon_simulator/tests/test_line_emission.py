from photon_simulator import \
    LineEmissionSourceModel, PhotonList
from photon_simulator.tests.beta_model_source import \
    BetaModelSource
from yt.units.yt_array import YTQuantity, uconcatenate
import numpy as np
import yt.units as u
from yt.utilities.physical_constants import clight

cross_section = 500.0e-22*u.cm**3/u.s
m_chi = (10.0*u.GeV).to_equivalent("g", "mass_energy")

def test_line_emission():

    bms = BetaModelSource()
    ds = bms.ds

    def _dm_emission(field, data):
        return cross_section*(data["dark_matter_density"]/m_chi)**2
    ds.add_field(("gas","dm_emission"), function=_dm_emission, units="s**-1*cm**-3")

    location = YTQuantity(3.5, "keV")
    sigma = YTQuantity(1000., "km/s")
    sigma_E = (location*sigma/clight).in_units("keV")

    A = YTQuantity(1000., "cm**2")
    exp_time = YTQuantity(2.0e5, "s")
    redshift = 0.01

    sphere = ds.sphere("c", (100.,"kpc"))

    line_model = LineEmissionSourceModel(location, "dm_emission", sigma="dark_matter_dispersion", prng=bms.prng)

    photons = PhotonList.from_data_source(sphere, redshift, A, exp_time,
                                          line_model)

    D_A = photons.parameters["FiducialAngularDiameterDistance"]
    dist_fac = 1.0/(4.*np.pi*D_A*D_A*(1.+redshift)**2)
    dm_E = (sphere["dm_emission"]*sphere["cell_volume"]).sum()

    E = uconcatenate(photons["Energy"])
    n_E = len(E)

    n_E_pred = (exp_time*A*dm_E*dist_fac).in_units("dimensionless")

    loc = location/(1.+redshift)
    sig = sigma_E/(1.+redshift)

    assert np.abs(loc-E.mean()) < sig/np.sqrt(n_E)
    assert np.abs(E.std()**2-sig*sig) < np.sqrt((3*(n_E-1)**2-(n_E-1)*(n_E-3))/n_E**3)*sig**2
    assert np.abs(n_E-n_E_pred) < np.sqrt(n_E)