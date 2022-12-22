"""
A unit test for the pyxsim analysis module.
"""

import os
import shutil
import tempfile

import numpy as np
from numpy.testing import assert_allclose
from soxs import ApecGenerator
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

    print(pvalue)
    assert pvalue > 0.05

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_beta_model_fields():
    bms = BetaModelSource()
    ds = bms.ds

    redshift = 0.2

    sphere = ds.sphere("c", (0.5, "Mpc"))

    kT_sim = bms.kT
    Z_sim = bms.Z

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")
    D_L = cosmo.luminosity_distance(0.0, redshift).to_value("cm")

    norm = (
        1.0e-14
        * sphere.sum(("gas", "emission_measure")).v
        / (4.0 * np.pi * D_A * D_A * (1 + redshift) ** 2)
    )

    agen = ApecGenerator(0.1, 11.5, 2000)

    spec = agen.get_spectrum(kT_sim, Z_sim, redshift, norm)
    pflux, eflux = spec.get_flux_in_band(0.5 / (1.0 + redshift), 7.0 / (1.0 + redshift))
    lum = 4.0 * np.pi * D_L**2 * eflux.value
    plum = 4.0 * np.pi * D_L**2 * pflux.value / (1.0 + redshift)

    thermal_model = CIESourceModel("apec", 0.1, 11.5, 2000, Z_sim)

    xray_fields = thermal_model.make_source_fields(ds, 0.5, 7.0)
    lum1 = sphere.sum(xray_fields[1]).value
    plum1 = (sphere[xray_fields[-1]] * sphere["cell_volume"]).sum().value

    int_fields = thermal_model.make_intensity_fields(
        ds, 0.5 / (1.0 + redshift), 7.0 / (1.0 + redshift), redshift=redshift
    )
    angular_scale = 1.0 / cosmo.angular_scale(0.0, redshift).to("cm/arcsec")

    eflux2 = (sphere[int_fields[0]] * sphere["cell_volume"]).sum() * angular_scale**2
    pflux2 = (sphere[int_fields[1]] * sphere["cell_volume"]).sum() * angular_scale**2

    assert np.abs(lum1 - lum) / lum < 0.001
    assert np.abs(plum1 - plum) / plum < 0.01

    assert np.abs(eflux2.value - eflux.value) / eflux.value < 0.001
    assert np.abs(pflux2.value - pflux.value) / pflux.value < 0.01


def test_beta_model_spectrum():
    bms = BetaModelSource()
    ds = bms.ds

    redshift = 0.2

    sphere = ds.sphere("c", (0.5, "Mpc"))

    kT_sim = bms.kT
    Z_sim = bms.Z

    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")

    norm1 = 1.0e-14 * sphere.sum(("gas", "emission_measure")).v
    norm2 = norm1 / (4.0 * np.pi * D_A * D_A * (1 + redshift) ** 2)

    agen = ApecGenerator(0.2, 7.0, 2000)

    spec1 = agen.get_spectrum(kT_sim, Z_sim, redshift, norm2)

    thermal_model = CIESourceModel("apec", 0.2, 7.0, 2000, Z_sim)
    spec2 = thermal_model.make_spectrum(
        sphere, 0.2, 7.0, 2000, redshift=redshift, cosmology=cosmo
    )
    assert_allclose(spec1.flux.value, spec2.flux.value)

    spec3 = agen.get_spectrum(kT_sim, Z_sim, 0.0, norm1)

    spec4 = thermal_model.make_spectrum(sphere, 0.2, 7.0, 2000)

    assert_allclose(spec3.flux.value, spec4.flux.value)
