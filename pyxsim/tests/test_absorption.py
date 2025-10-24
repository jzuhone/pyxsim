import os
import shutil
import tempfile
from copy import deepcopy

import h5py
import numpy as np
import soxs
from yt.utilities.cosmology import Cosmology

from pyxsim import (
    CIESourceModel,
    make_column_density_map,
    make_photons,
    project_photons,
)
from pyxsim.tests.utils import UniformSource

cosmo = Cosmology()


def test_absorption():
    us = UniformSource()
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = us.prng

    emin = 0.1
    emax = 11.5
    nspec = 10000
    A = 3000.0
    exp_time = 1.0e4
    redshift = 0.05
    nH_sim = 0.02

    cosmo = Cosmology()
    D_A = cosmo.angular_diameter_distance(0.0, redshift).to_value("cm")

    dd = us.ds.all_data()

    thermal_model = CIESourceModel("apec", emin, emax, nspec, us.Z, prng=prng)
    make_photons("my_photons", dd, redshift, A, exp_time, thermal_model)
    make_column_density_map(
        us.ds,
        "z",
        "c",
        (1.0, "unitary"),
        (1.0, "unitary"),
        us.nx,
        us.nx,
        "nH.h5",
        field=("gas", "H_p0_number_density"),
    )
    project_photons(
        "my_photons",
        "my_events",
        "z",
        [30.0, 45.0],
        absorb_model="tbabs",
        nH=nH_sim,
        prng=prng,
        column_file="nH.h5",
        save_los=True,
    )
    with h5py.File("nH.h5", "r") as f:
        nH_cube = {key: f["data"][key][:] for key in f["data"].keys()}
    with h5py.File("my_events.h5", "r") as f:
        events = {key: f["data"][key][:] for key in f["data"].keys()}
    agen = soxs.ApecGenerator(emin, emax, nspec)
    dd = us.ds.all_data()
    EM = dd.sum(("gas", "emission_measure")) / us.nx
    norm = (
        exp_time
        * A
        * 1.0e-14
        * EM.in_base("cgs")
        / (4.0 * np.pi * D_A * D_A * (1.0 + redshift) ** 2)
    )
    specm = agen.get_spectrum(us.kT, us.Z, redshift, norm)
    specm.apply_foreground_absorption(nH_sim, model="tbabs")
    for i in range(0, us.nx - 1, 10):
        wbin = (nH_cube["wbins"][i], nH_cube["wbins"][i + 1])
        these_events = (events["los"] >= wbin[0]) & (events["los"] < wbin[1])
        eobs = events["eobs"][these_events]
        spec, _ = np.histogram(eobs, bins=specm.ebins.value)
        specm_i = deepcopy(specm)
        specm_i.apply_foreground_absorption(
            nH_cube["nH"][0, 0, i], model="tbabs", redshift=redshift
        )
        with np.errstate(divide="ignore"):
            err = (specm_i.binned_flux.value - spec) / np.sqrt(spec)
        err[spec == 0.0] = 0.0
        chi2 = (err**2).sum() / nspec
        assert chi2 < 1.2
    os.chdir(curdir)
    shutil.rmtree(tmpdir)
