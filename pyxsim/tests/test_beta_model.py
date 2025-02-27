"""
A unit test for the pyxsim analysis module.
"""

import os
import shutil
import tempfile

import numpy as np
from numpy.testing import assert_allclose
from soxs import ApecGenerator
from yt import YTQuantity
from yt.utilities.cosmology import Cosmology
from yt.utilities.physical_constants import clight

from pyxsim import CIESourceModel, make_photons, project_photons
from pyxsim.tests.utils import (
    BetaModelSource,
    ParticleBetaModelSource,
    events_ks_testing,
)

cosmo = Cosmology()

ckms = clight.in_units("km/s").v


def test_beta_model(check_dir):
    bms = BetaModelSource()
    do_beta_model(bms, check_dir)


def test_beta_model_nomove(check_dir):
    bms = BetaModelSource()
    do_beta_model(bms, check_dir, axis="x", prng=89)


def test_beta_model_offaxis(check_dir):
    bms = BetaModelSource()
    do_beta_model(bms, check_dir, axis=[1.0, -2.0, 5.0], prng=78)


def test_particle_beta_model(check_dir):
    bms = ParticleBetaModelSource()
    do_beta_model(bms, check_dir, prng=29)


def test_particle_beta_model_nomove(check_dir):
    bms = ParticleBetaModelSource()
    do_beta_model(bms, check_dir, axis="x", prng=72)


def test_particle_beta_model_offaxis(check_dir):
    bms = ParticleBetaModelSource()
    do_beta_model(bms, check_dir, prng=67, axis=[1.0, -2.0, 5.0])


def do_beta_model(source, check_dir, axis="z", prng=None):
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    if prng is None:
        prng = source.prng

    ds = source.ds

    A = 30000.0
    exp_time = 1.0e4
    redshift = 0.05
    nH_sim = 0.02

    sphere = ds.sphere("c", (0.5, "Mpc"))

    kT_sim = source.kT
    Z_sim = source.Z

    thermal_model = CIESourceModel("apec", 0.1, 11.5, 20000, Z_sim, prng=prng)
    make_photons("my_photons", sphere, redshift, A, exp_time, thermal_model)

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")

    norm_sim = sphere.quantities.total_quantity(("gas", "emission_measure"))
    norm_sim *= 1.0e-14 / (4 * np.pi * D_A * D_A * (1.0 + redshift) * (1.0 + redshift))
    norm_sim = float(norm_sim.in_cgs())

    v1, v2 = sphere.quantities.weighted_standard_deviation(
        ("gas", "velocity_z"), ("gas", "emission_measure")
    )

    if isinstance(axis, str):
        if axis == "z":
            fac = 1.0
        else:
            fac = 0.0
    else:
        axis /= np.sqrt(np.dot(axis, axis))
        fac = np.dot(axis, [0.0, 0.0, 1.0])

    sigma_sim = fac * float(v1.in_units("km/s"))
    mu_sim = -fac * float(v2.in_units("km/s"))

    project_photons(
        "my_photons",
        "my_events",
        axis,
        [30.0, 45.0],
        absorb_model="tbabs",
        nH=nH_sim,
        prng=prng,
    )

    redshift_sim = (1.0 + mu_sim / ckms) * (1.0 + redshift) - 1.0

    agen = ApecGenerator(0.3, 8.0, 10000)
    spec = agen.get_spectrum(kT_sim, Z_sim, redshift_sim, norm_sim, velocity=sigma_sim)
    spec.apply_foreground_absorption(nH_sim, model="tbabs")

    pvalue = events_ks_testing("my_events.h5", spec, exp_time, A, check_dir)

    assert pvalue > 0.05

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_vapec_beta_model(check_dir):
    bms = BetaModelSource()

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = 50

    ds = bms.ds

    A = 30000.0
    exp_time = 1.0e5
    redshift = 0.05
    nH_sim = 0.02

    sphere = ds.sphere("c", (0.5, "Mpc"))

    kT_sim = bms.kT
    Z_sim = bms.Z
    O_sim = bms.O
    Ca_sim = bms.Ca

    var_elem = {"O": ("stream", "oxygen"), "Ca": ("stream", "calcium")}

    thermal_model = CIESourceModel(
        "apec", 0.1, 11.5, 20000, ("gas", "metallicity"), var_elem=var_elem, prng=prng
    )

    make_photons("my_photons", sphere, redshift, A, exp_time, thermal_model)

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to("cm")

    norm_sim = sphere.quantities.total_quantity("emission_measure")
    norm_sim *= 1.0e-14 / (4 * np.pi * D_A * D_A * (1.0 + redshift) * (1.0 + redshift))
    norm_sim = float(norm_sim.in_cgs())

    project_photons(
        "my_photons",
        "my_events",
        "z",
        [30.0, 45.0],
        absorb_model="tbabs",
        nH=nH_sim,
        prng=prng,
        no_shifting=True,
    )

    agen = ApecGenerator(0.2, 10.0, 10000, var_elem=["O", "Ca"])
    spec = agen.get_spectrum(
        kT_sim, Z_sim, redshift, norm_sim, elem_abund={"O": O_sim, "Ca": Ca_sim}
    )
    spec.apply_foreground_absorption(nH_sim, model="tbabs")

    pvalue = events_ks_testing("my_events.h5", spec, exp_time, A, check_dir)

    assert pvalue > 0.05

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_beta_model_fields():
    from astropy.cosmology import FlatLambdaCDM

    vlos = -0.5
    vtot = np.abs(vlos)
    v_s = YTQuantity(vlos, "c").to_value("cm/s")
    bms = BetaModelSource(no_broad=True, v_s=v_s)
    ds = bms.ds

    shift = np.sqrt(1.0 - vtot**2) / (1.0 - vlos)

    redshift = 0.2

    acosmo = FlatLambdaCDM(H0=71.0, Om0=0.27)

    kT_sim = bms.kT
    Z_sim = bms.Z

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")

    sphere = ds.sphere("c", (0.2, "Mpc"))

    norm = (
        1.0e-14
        * sphere.sum(("gas", "emission_measure")).v
        / (4.0 * np.pi * D_A * D_A * (1 + redshift) ** 2)
    )

    agen = ApecGenerator(0.1, 11.5, 10000)

    spec = agen.get_spectrum(kT_sim, Z_sim, redshift, norm)
    pflux, eflux = spec.get_flux_in_band(1.0, 4.0)
    plum, lum = spec.get_lum_in_band(1.0, 4.0, redshift=redshift, cosmology=acosmo)

    thermal_model = CIESourceModel("apec", 0.1, 11.5, 10000, Z_sim)

    xray_fields = thermal_model.make_source_fields(ds, 1.0, 4.0)
    lum1 = (sphere[xray_fields[0]] * sphere["cell_volume"]).sum().value
    lum2 = sphere.sum(xray_fields[1]).value
    plum1 = (sphere[xray_fields[-2]] * sphere["cell_volume"]).sum().value
    plum2 = sphere[xray_fields[-1]].sum().value

    assert np.abs(lum1 - lum.value) / lum.value < 0.001
    assert np.abs(lum2 - lum.value) / lum.value < 0.001
    assert np.abs(plum1 - plum.value) / plum.value < 0.0012
    assert np.abs(plum2 - plum.value) / plum.value < 0.0012

    angular_scale = 1.0 / cosmo.angular_scale(0.0, redshift).to("cm/arcsec")

    sphere2 = ds.sphere("c", (0.2, "Mpc"))

    sphere2.set_field_parameter("axis", 2)

    int_fields = thermal_model.make_intensity_fields(
        ds,
        1.0,
        4.0,
        redshift=redshift,
    )

    eflux3 = (sphere2[int_fields[0]] * sphere2["cell_volume"]).sum() * angular_scale**2
    pflux3 = (sphere2[int_fields[1]] * sphere2["cell_volume"]).sum() * angular_scale**2

    pflux4, eflux4 = spec.get_flux_in_band(1.0 / shift, 4.0 / shift)

    eflux4 *= shift**4
    pflux4 *= shift**3

    assert np.abs(eflux3.value - eflux4.value) / eflux4.value < 0.001
    assert np.abs(pflux3.value - pflux4.value) / pflux4.value < 0.001

    sphere3 = ds.sphere("c", (0.2, "Mpc"))

    sphere3.set_field_parameter("axis", 2)

    int_fields = thermal_model.make_intensity_fields(
        ds,
        1.0,
        4.0,
        redshift=redshift,
        no_doppler=True,
        force_override=True,
        band_name="no_shift",
    )

    eflux5 = (sphere3[int_fields[0]] * sphere3["cell_volume"]).sum() * angular_scale**2
    pflux5 = (sphere3[int_fields[1]] * sphere3["cell_volume"]).sum() * angular_scale**2

    assert np.abs(eflux5.value - eflux.value) / eflux.value < 0.001
    assert np.abs(pflux5.value - pflux.value) / pflux.value < 0.001


def test_beta_model_spectrum():

    vlos = -0.2
    v_s = YTQuantity(vlos, "c").to_value("cm/s")
    bms = BetaModelSource(no_broad=True, v_s=v_s)
    ds = bms.ds

    redshift = 0.2

    sphere = ds.sphere("c", (0.5, "Mpc"))

    kT_sim = bms.kT
    Z_sim = bms.Z

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")

    norm1 = 1.0e-14 * sphere.sum(("gas", "emission_measure")).v
    norm2 = norm1 / (4.0 * np.pi * D_A * D_A * (1 + redshift) ** 2)

    agen = ApecGenerator(0.1, 11.0, 3000)

    spec1 = agen.get_spectrum(kT_sim, Z_sim, redshift, norm2).regrid_spectrum(
        0.2, 6.0, 2000
    )

    thermal_model = CIESourceModel("apec", 0.1, 11.0, 3000, Z_sim)
    spec2 = thermal_model.make_spectrum(
        sphere, 0.2, 6.0, 2000, redshift=redshift, cosmology=cosmo
    )
    assert_allclose(spec1.flux.value, spec2.flux.value)

    spec3 = agen.get_spectrum(kT_sim, Z_sim, 0.0, norm1).regrid_spectrum(0.2, 6.0, 2000)

    spec4 = thermal_model.make_spectrum(sphere, 0.2, 6.0, 2000)

    assert_allclose(spec3.flux.value, spec4.flux.value)

    spec5 = agen.get_spectrum(kT_sim, Z_sim, redshift, norm2).regrid_spectrum(
        0.2, 6.0, 2000, vlos=-vlos * ckms
    )

    spec6 = thermal_model.make_spectrum(
        sphere,
        0.2,
        6.0,
        2000,
        redshift=redshift,
        cosmology=cosmo,
        normal="z",
    )
    assert_allclose(spec5.flux.value, spec6.flux.value)

    spec7 = agen.get_spectrum(kT_sim, Z_sim, redshift, norm2).regrid_spectrum(
        0.2, 6.0, 2000, vlos=0.0, vtot=vlos * ckms
    )

    spec8 = thermal_model.make_spectrum(
        sphere,
        0.2,
        6.0,
        2000,
        redshift=redshift,
        cosmology=cosmo,
        normal="x",
    )
    assert_allclose(spec7.flux.value, spec8.flux.value)
