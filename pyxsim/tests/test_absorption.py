import os
import shutil
import tempfile

from yt.utilities.cosmology import Cosmology

from pyxsim import (
    CIESourceModel,
    make_column_density_map,
    make_photons,
    project_photons,
)
from pyxsim.tests.utils import UniformSource

cosmo = Cosmology()


def test_absorption(prng=None):
    us = UniformSource()
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    if prng is None:
        prng = us.prng

    A = 10000.0
    exp_time = 1.0e4
    redshift = 0.05
    nH_sim = 0.02

    dd = us.ds.all_data()

    thermal_model = CIESourceModel("apec", 0.1, 11.5, 20000, us.Z, prng=prng)
    make_photons("my_photons", dd, redshift, A, exp_time, thermal_model)
    make_column_density_map(
        us.ds,
        "z",
        "c",
        (1.0, "unitary"),
        (1.0, "unitary"),
        100,
        100,
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
    )
    os.system("cp nH.h5 /Users/jzuhone")
    os.chdir(curdir)
    shutil.rmtree(tmpdir)


if __name__ == "__main__":
    test_absorption()
