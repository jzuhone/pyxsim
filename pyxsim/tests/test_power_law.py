import os
import shutil
import tempfile

import numpy as np
from numpy.testing import assert_allclose
from soxs import CountRateSpectrum, Spectrum
from soxs.constants import keV_per_erg
from yt.units.yt_array import YTQuantity
from yt.utilities.cosmology import Cosmology
from yt.utilities.physical_constants import mp

from pyxsim import PowerLawSourceModel, make_photons, project_photons
from pyxsim.tests.utils import BetaModelSource, events_ks_testing


def test_power_law(check_dir):
    plaw_fit(1.1, check_dir, prng=33)
    plaw_fit(0.8, check_dir, prng=28)
    plaw_fit(1.0, check_dir, prng=23)


def plaw_fit(alpha_sim, check_dir, prng=None):
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    bms = BetaModelSource()
    ds = bms.ds

    if prng is None:
        prng = bms.prng

    def _power_law_luminosity(field, data):
        norm = data.ds.quan(
            1.3e-26, "erg/s"
        )  # this seems small, but only because it is per-cell
        return norm * data["cell_mass"] / (1.0 * mp)

    ds.add_field(
        ("gas", "plaw_lum"),
        function=_power_law_luminosity,
        sampling_type="local",
        units="erg/s",
    )

    nH_sim = 0.02

    A = YTQuantity(2000.0, "cm**2")
    exp_time = YTQuantity(2.0e5, "s")
    redshift = 0.01

    cosmo = Cosmology()

    sphere = ds.sphere("c", (100.0, "kpc"))

    plaw_model = PowerLawSourceModel(
        1.0, 0.01, 11.0, ("gas", "plaw_lum"), alpha_sim, prng=prng
    )

    make_photons("plaw_photons.h5", sphere, redshift, A, exp_time, plaw_model)

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")
    dist_fac = 1.0 / (4.0 * np.pi * D_A * D_A * (1.0 + redshift) ** 3)
    norm_sim = float((sphere["plaw_lum"].to("keV/s")).sum() * dist_fac) * (
        1.0 + redshift
    )
    norm_sim *= (2.0 - alpha_sim) / (
        11.0 ** (2.0 - alpha_sim) - 0.01 ** (2.0 - alpha_sim)
    )

    project_photons(
        "plaw_photons.h5",
        "plaw_events.h5",
        "z",
        [30.0, 45.0],
        absorb_model="wabs",
        nH=nH_sim,
        prng=bms.prng,
        no_shifting=True,
    )

    spec = Spectrum.from_powerlaw(alpha_sim, redshift, norm_sim, 0.2, 10.0, 5000)
    spec.apply_foreground_absorption(nH_sim, model="wabs")

    pvalue = events_ks_testing(
        "plaw_events.h5", spec, exp_time.value, A.value, check_dir
    )

    assert pvalue > 0.05

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_power_law_fields():
    bms = BetaModelSource()
    ds = bms.ds

    def _power_law_luminosity(field, data):
        norm = data.ds.quan(
            1.3e-26, "erg/s"
        )  # this seems small, but only because it is per-cell
        return norm * data["cell_mass"] / (1.0 * mp)

    ds.add_field(
        ("gas", "hard_emission"),
        function=_power_law_luminosity,
        sampling_type="local",
        units="erg/s",
    )

    sphere = ds.sphere("c", (100.0, "kpc"))

    lum = sphere["gas", "hard_emission"].sum().v

    alpha1 = 1.1
    plaw_model1 = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", alpha1)

    src_fields1 = plaw_model1.make_source_fields(ds, 0.5, 7.0)

    elum0 = (
        lum
        * (7.0 ** (2.0 - alpha1) - 0.5 ** (2.0 - alpha1))
        / (11.0 ** (2.0 - alpha1) - 0.01 ** (2.0 - alpha1))
    )
    elum1 = sphere[src_fields1[1]].sum().to_value("erg/s")
    assert_allclose(elum0, elum1)

    plum0 = (
        keV_per_erg
        * lum
        * (7.0 ** (1.0 - alpha1) - 0.5 ** (1.0 - alpha1))
        / (11.0 ** (2.0 - alpha1) - 0.01 ** (2.0 - alpha1))
        * (2.0 - alpha1)
        / (1.0 - alpha1)
    )
    plum1 = (sphere[src_fields1[-2]] * sphere["cell_volume"]).sum().v
    plumr = sphere[src_fields1[-1]].sum().to_value("photons/s")
    assert_allclose(plum0, plum1)
    assert_allclose(plum0, plumr)

    del sphere[src_fields1[-1]]
    del sphere[src_fields1[-2]]
    del sphere[src_fields1[1]]

    alpha2 = 1.0
    plaw_model2 = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", alpha2)

    src_fields2 = plaw_model2.make_source_fields(ds, 0.5, 7.0, force_override=True)

    elum2 = lum * (7.0 - 0.5) / (11.0 - 0.01)
    elum3 = sphere[src_fields2[1]].sum().to_value("erg/s")
    assert_allclose(elum2, elum3)

    plum2 = (
        keV_per_erg
        * lum
        * np.log(7.0 / 0.5)
        / (11.0 ** (2.0 - alpha2) - 0.01 ** (2.0 - alpha2))
        * (2.0 - alpha2)
    )
    plum3 = (sphere[src_fields2[-2]] * sphere["cell_volume"]).sum().v
    plum4 = sphere[src_fields2[-1]].sum().to_value("photons/s")
    assert_allclose(plum2, plum3)
    assert_allclose(plum2, plum4)

    del sphere[src_fields2[-1]]
    del sphere[src_fields2[-2]]
    del sphere[src_fields2[1]]

    alpha3 = 2.0
    plaw_model3 = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", alpha3)

    src_fields3 = plaw_model3.make_source_fields(ds, 0.5, 7.0, force_override=True)

    elum4 = lum * np.log(7.0 / 0.5) / np.log(11.0 / 0.01)
    elum5 = sphere[src_fields3[1]].sum().to_value("erg/s")
    assert_allclose(elum4, elum5)

    plum5 = (
        keV_per_erg
        * lum
        * (7.0 ** (1.0 - alpha3) - 0.5 ** (1.0 - alpha3))
        / np.log(11.0 / 0.01)
        / (1.0 - alpha3)
    )
    plum6 = (sphere[src_fields3[-2]] * sphere["cell_volume"]).sum().v
    plum7 = sphere[src_fields3[-1]].sum().to_value("photons/s")
    assert_allclose(plum5, plum6)
    assert_allclose(plum5, plum7)


def test_power_law_spectrum():
    cosmo = Cosmology()

    vtot = vlos = 0.5
    v_s = YTQuantity(vlos, "c").to_value("cm/s")
    bms = BetaModelSource(no_broad=True, v_s=v_s)

    ds = bms.ds

    redshift = 0.2

    def _power_law_luminosity(field, data):
        norm = data.ds.quan(
            1.3e-26, "erg/s"
        )  # this seems small, but only because it is per-cell
        return norm * data["cell_mass"] / (1.0 * mp)

    ds.add_field(
        ("gas", "hard_emission"),
        function=_power_law_luminosity,
        sampling_type="local",
        units="erg/s",
    )

    sphere = ds.sphere("c", (100.0, "kpc"))
    lum = sphere["gas", "hard_emission"].sum().v

    alpha = 1.1
    plaw_model = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", alpha)

    lum_scaled = (
        lum
        * (6.0 ** (2.0 - alpha) - 0.1 ** (2.0 - alpha))
        / (11.0 ** (2.0 - alpha) - 0.01 ** (2.0 - alpha))
    )
    spec0 = CountRateSpectrum.from_powerlaw(alpha, 0.0, 1.0, 0.1, 6.0, 1000)
    spec0.rescale_flux(lum_scaled, flux_type="energy")
    spec1 = plaw_model.make_spectrum(sphere, 0.1, 6.0, 1000)

    assert_allclose(spec0.flux.value, spec1.flux.value, rtol=1.0e-6)

    D_L = cosmo.luminosity_distance(0.0, redshift).to_value("cm")

    lum_scaled2 = (
        lum
        * (
            (6.0 * (1.0 + redshift)) ** (2.0 - alpha)
            - (0.1 * (1.0 + redshift)) ** (2.0 - alpha)
        )
        / (11.0 ** (2.0 - alpha) - 0.01 ** (2.0 - alpha))
    )
    spec2 = Spectrum.from_powerlaw(alpha, redshift, 1.0, 0.1, 6.0, 1000)
    spec2.rescale_flux(lum_scaled2 / (4.0 * np.pi * D_L**2), flux_type="energy")
    spec3 = plaw_model.make_spectrum(sphere, 0.1, 6.0, 1000, redshift=redshift)

    assert_allclose(spec2.flux.value, spec3.flux.value, rtol=1.0e-6)

    shift = np.sqrt(1.0 - vtot**2) / (1.0 - vlos)
    spec4 = plaw_model.make_spectrum(
        sphere, 0.1, 6.0, 1000, redshift=redshift, normal=[0.0, 0.0, 1.0]
    )

    assert_allclose(
        spec2.flux.value * (shift ** (2.0 + alpha)), spec4.flux.value, rtol=1.0e-6
    )
