import os
import shutil
import tempfile

import h5py
import numpy as np
from numpy.testing import assert_allclose
from unyt import clight
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.cosmology import Cosmology

from pyxsim import LineSourceModel, make_photons
from pyxsim.tests.utils import BetaModelSource

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
        return (
            (data["gas", "dark_matter_density"] / m_chi) ** 2
            * data["cell_volume"]
            * cross_section
        )

    ds.add_field(
        ("gas", "dm_emission"),
        function=_dm_emission,
        units="s**-1",
        sampling_type="cell",
    )

    location = YTQuantity(3.5, "keV")
    sigma = YTQuantity(1000.0, "km/s")
    sigma_E = (location * sigma / clight).in_units("keV")

    A = YTQuantity(1000.0, "cm**2")
    exp_time = YTQuantity(2.0e5, "s")
    redshift = 0.01

    sphere = ds.sphere("c", (100.0, "kpc"))

    line_model = LineSourceModel(
        location, "dm_emission", sigma=("stream", "dark_matter_dispersion"), prng=32
    )

    n_photons, n_cells = make_photons(
        "my_photons.h5", sphere, redshift, A, exp_time, line_model
    )

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to("cm")

    dist_fac = 1.0 / (4.0 * np.pi * D_A * D_A * (1.0 + redshift) ** 3)
    dm_E = (sphere["dm_emission"]).sum()

    with h5py.File("my_photons.h5", "r") as f:
        E = YTArray(f["data"]["energy"][:], "keV")
        n_E = len(E)

    n_E_pred = (exp_time * A * dm_E * dist_fac).in_units("dimensionless")

    loc = location / (1.0 + redshift)
    sig = sigma_E / (1.0 + redshift)

    assert np.abs(loc - E.mean()) < 1.645 * sig / np.sqrt(n_E)
    assert (
        np.abs(E.std() ** 2 - sig * sig) < 1.645 * np.sqrt(2 * (n_E - 1)) * sig**2 / n_E
    )
    assert np.abs(n_E - n_E_pred) < 1.645 * np.sqrt(n_E)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_line_emission_fields():
    cosmo = Cosmology()

    vtot = vlos = 0.5
    v_s = YTQuantity(vlos, "c").to_value("cm/s")
    bms = BetaModelSource(no_broad=True, v_s=v_s)
    ds = bms.ds

    redshift = 0.2

    angular_scale = 1.0 / cosmo.angular_scale(0.0, redshift).to("cm/arcsec")

    def _dm_emission(field, data):
        return (
            (data["gas", "dark_matter_density"] / m_chi) ** 2
            * data["cell_volume"]
            * cross_section
        )

    ds.add_field(
        ("gas", "dm_emission"),
        function=_dm_emission,
        units="s**-1",
        sampling_type="cell",
    )

    location = YTQuantity(3.5, "keV")
    sigma = YTQuantity(1000.0, "km/s")
    sigma_E = (location * sigma / clight).in_units("keV")

    sphere = ds.sphere("c", (100.0, "kpc"))

    dm_E = (sphere["dm_emission"]).sum()

    line_model1 = LineSourceModel(location, "dm_emission", (1.0, "eV"))

    line_fields1 = line_model1.make_source_fields(ds, 0.5, 7.0)
    assert_allclose(sphere[line_fields1[1]].sum().to("keV/s"), dm_E * location)
    assert_allclose((sphere[line_fields1[-2]] * sphere["cell_volume"]).sum(), dm_E)
    assert_allclose(sphere[line_fields1[-1]].sum(), dm_E)

    line_fields2 = line_model1.make_source_fields(ds, 0.5, 2.0)
    for field in line_fields2:
        assert_allclose(sphere[field].sum(), 0.0)

    line_model2 = LineSourceModel(
        location, "dm_emission", ("stream", "dark_matter_dispersion")
    )

    del sphere[line_fields1[-2]]
    del sphere[line_fields1[1]]

    line_fields3 = line_model2.make_source_fields(ds, 0.5, 7.0, force_override=True)
    assert_allclose(sphere[line_fields3[1]].sum().to("keV/s"), dm_E * location)
    assert_allclose((sphere[line_fields3[-2]] * sphere["cell_volume"]).sum(), dm_E)

    line_fields4 = line_model2.make_source_fields(ds, 0.1, 3.5)
    de = sigma_E / np.sqrt(2.0 * np.pi)

    assert_allclose(
        sphere[line_fields4[1]].sum().to("keV/s"), dm_E * (0.5 * location - de)
    )
    assert_allclose(
        (sphere[line_fields4[-2]] * sphere["cell_volume"]).sum(), 0.5 * dm_E
    )
    assert_allclose(sphere[line_fields4[-1]].sum(), 0.5 * dm_E)

    sphere.set_field_parameter("axis", 2)

    int_fields = line_model2.make_intensity_fields(ds, 0.5, 7.0, redshift=redshift)

    dist_fac = (
        1.0 / (4.0 * np.pi * cosmo.luminosity_distance(0.0, redshift).to("cm") ** 2).v
    )
    shift = np.sqrt(1.0 - vtot**2) / (1.0 - vlos)

    eflux = (sphere[int_fields[0]] * sphere["cell_volume"]).sum() * angular_scale**2
    pflux = (sphere[int_fields[1]] * sphere["cell_volume"]).sum() * angular_scale**2

    assert_allclose(eflux.v, (shift**4) * location.to_value("erg") * dm_E.v * dist_fac)
    assert_allclose(pflux.v, (1.0 + redshift) * (shift**3) * dm_E.v * dist_fac)


def test_line_emission_spectra():
    cosmo = Cosmology()

    vtot = vlos = 0.5
    v_s = YTQuantity(vlos, "c").to_value("cm/s")
    bms = BetaModelSource(no_broad=True, v_s=v_s)
    ds = bms.ds

    redshift = 0.2

    def _dm_emission(field, data):
        return data["index", "ones"] * data.ds.quan(1.0e40, "1/s")

    ds.add_field(
        ("gas", "dm_emission"),
        function=_dm_emission,
        units="s**-1",
        sampling_type="cell",
    )

    sphere = ds.sphere("c", (100.0, "kpc"))

    dm_E = (sphere["dm_emission"]).sum()

    location = YTQuantity(3.5, "keV")
    sigma = YTQuantity(1000.0, "km/s")
    sigma_E = (location * sigma / clight).in_units("keV")

    line_model = LineSourceModel(location, "dm_emission", sigma)

    def weight_std(x, w):
        return np.sqrt(np.average((x - np.average(x, weights=w)) ** 2, weights=w))

    spec1 = line_model.make_spectrum(sphere, 1.0, 6.0, 5000)
    assert_allclose(np.average(spec1.emid.value, weights=spec1.flux.value), location.v)
    assert_allclose(
        weight_std(spec1.emid.value, spec1.flux.value), sigma_E.v, rtol=1.0e-3
    )
    assert_allclose(np.sum(spec1.flux.value * spec1.de.value), dm_E.v)

    spec2 = line_model.make_spectrum(sphere, 1.0, 6.0, 5000, redshift=redshift)
    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")

    dist_fac = 1.0 / (4.0 * np.pi * (D_A * (1.0 + redshift)) ** 2)

    assert_allclose(
        np.average(spec2.emid.value, weights=spec2.flux.value),
        location.v / (1.0 + redshift),
    )
    assert_allclose(
        weight_std(spec2.emid.value, spec2.flux.value),
        sigma_E.v / (1.0 + redshift),
        rtol=1.0e-4,
    )
    assert_allclose(
        np.sum(spec2.flux.value * spec2.de.value), dm_E.v * dist_fac / (1.0 + redshift)
    )

    spec3 = line_model.make_spectrum(
        sphere, 1.0, 6.0, 5000, redshift=redshift, normal=[0.0, 0.0, 1.0]
    )

    shift = np.sqrt(1.0 - vtot**2) / (1.0 - vlos)

    assert_allclose(
        np.average(spec3.emid.value, weights=spec3.flux.value),
        location.v * shift / (1.0 + redshift),
    )

    assert_allclose(
        weight_std(spec3.emid.value, spec3.flux.value),
        sigma_E.v * shift / (1.0 + redshift),
        rtol=1.0e-4,
    )

    assert_allclose(
        np.sum(spec3.flux.value * spec3.de.value),
        dm_E.v * dist_fac * shift**3 / (1.0 + redshift),
    )
