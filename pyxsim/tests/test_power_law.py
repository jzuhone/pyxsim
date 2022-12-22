import os
import shutil
import tempfile

import numpy as np
from numpy.testing import assert_allclose
from soxs import Spectrum
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

    def _hard_emission(field, data):
        return (
            data.ds.quan(1.0e-18, "s**-1*keV**-1")
            * data["density"]
            * data["cell_volume"]
            / mp
        )

    ds.add_field(
        ("gas", "hard_emission"),
        function=_hard_emission,
        units="keV**-1*s**-1",
        sampling_type="local",
    )

    nH_sim = 0.02

    A = YTQuantity(2000.0, "cm**2")
    exp_time = YTQuantity(2.0e5, "s")
    redshift = 0.01

    cosmo = Cosmology()

    sphere = ds.sphere("c", (100.0, "kpc"))

    plaw_model = PowerLawSourceModel(
        1.0, 0.01, 11.0, "hard_emission", alpha_sim, prng=prng
    )

    make_photons("plaw_photons.h5", sphere, redshift, A, exp_time, plaw_model)

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")
    dist_fac = 1.0 / (4.0 * np.pi * D_A * D_A * (1.0 + redshift) ** 3)
    norm_sim = float((sphere["hard_emission"]).sum() * dist_fac) * (1.0 + redshift)

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

    def _hard_emission(field, data):
        return (
            data.ds.quan(1.0e-18, "s**-1*keV**-1")
            * data["density"]
            * data["cell_volume"]
            / mp
        )

    ds.add_field(
        ("gas", "hard_emission"),
        function=_hard_emission,
        units="keV**-1*s**-1",
        sampling_type="local",
    )

    sphere = ds.sphere("c", (100.0, "kpc"))

    norm = sphere["gas", "hard_emission"].sum().v

    alpha1 = 1.1
    plaw_model1 = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", alpha1)

    src_fields1 = plaw_model1.make_source_fields(ds, 0.5, 7.0)

    plum0 = -norm * (7.0**-0.1 - 0.5**-0.1) / 0.1
    plum1 = (sphere[src_fields1[-1]] * sphere["cell_volume"]).sum().v
    assert_allclose(plum0, plum1)

    elum0 = norm * (7.0**0.9 - 0.5**0.9) / 0.9
    elum1 = sphere[src_fields1[1]].sum().to_value("keV/s")
    assert_allclose(elum0, elum1)

    del sphere[src_fields1[-1]]
    del sphere[src_fields1[1]]

    alpha2 = 1.0
    plaw_model2 = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", alpha2)

    src_fields2 = plaw_model2.make_source_fields(ds, 0.5, 7.0, force_override=True)

    plum2 = norm * np.log(7.0 / 0.5)
    plum3 = (sphere[src_fields2[-1]] * sphere["cell_volume"]).sum().v
    assert_allclose(plum2, plum3)

    elum2 = norm * (7.0 - 0.5)
    elum3 = sphere[src_fields2[1]].sum().to_value("keV/s")
    assert_allclose(elum2, elum3)


def test_power_law_spectrum():

    cosmo = Cosmology()

    bms = BetaModelSource()
    ds = bms.ds

    def _hard_emission(field, data):
        return (
            data.ds.quan(1.0e-18, "s**-1*keV**-1")
            * data["density"]
            * data["cell_volume"]
            / mp
        )

    ds.add_field(
        ("gas", "hard_emission"),
        function=_hard_emission,
        units="keV**-1*s**-1",
        sampling_type="local",
    )

    sphere = ds.sphere("c", (100.0, "kpc"))
    norm = sphere["gas", "hard_emission"].sum().v

    alpha = 1.1
    plaw_model = PowerLawSourceModel(1.0, 0.01, 11.0, "hard_emission", alpha)

    spec0 = Spectrum.from_powerlaw(alpha, 0.0, norm, 0.1, 6.0, 1000)
    spec1 = plaw_model.make_spectrum(sphere, 0.1, 6.0, 1000)

    assert_allclose(spec0.flux.value, spec1.flux.value)

    D_A = cosmo.angular_diameter_distance(0.0, 0.2).to_value("cm")

    dist_fac = 1.0 / (4.0 * np.pi * (D_A * 1.2) ** 2)
    spec2 = Spectrum.from_powerlaw(alpha, 0.2, norm * dist_fac, 0.1, 6.0, 1000)
    spec3 = plaw_model.make_spectrum(sphere, 0.1, 6.0, 1000, redshift=0.2)

    assert_allclose(spec2.flux.value, spec3.flux.value)
