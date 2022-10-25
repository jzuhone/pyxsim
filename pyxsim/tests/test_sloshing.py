"""
Answer test pyxsim.
"""
import numpy as np

from pyxsim import \
    CIESourceModel, \
    EventList, make_photons, \
    project_photons, merge_files
from pyxsim.tests.utils import hdf5_answer_testing, file_answer_testing
from numpy.testing import assert_array_equal
from numpy.random import RandomState
import yt
import h5py

gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"


def test_sloshing(answer_store, answer_dir):

    prng = RandomState(0x4d3d3d3)

    ds = yt.load(gslr, default_species_fields="ionized")
    A = 2000.
    exp_time = 1.0e4
    redshift = 0.1

    sphere = ds.sphere("c", (0.5, "Mpc"))
    sphere.set_field_parameter("X_H", 0.75)

    thermal_model = CIESourceModel("apec", 0.1, 11.0, 10000, 0.3,
                                   thermal_broad=False, prng=prng)
    n_photons1, n_cells1 = make_photons("photons1", sphere, redshift, A, 
                                        exp_time, thermal_model)
 
    hdf5_answer_testing("photons1.h5", answer_store, answer_dir)

    n_events1 = project_photons("photons1", "events1", [1.0, -0.5, 0.2],
                                [30., 45.], absorb_model="tbabs", nH=0.1,
                                prng=prng)

    hdf5_answer_testing("events1.h5", answer_store, answer_dir)

    events1 = EventList("events1.h5")

    events1.write_fits_file("test_events.fits", (20.0, "arcmin"), 1024)
    events1.write_spectrum("test_spec.fits", 0.2, 10.0, 2000)
    events1.write_fits_image("test_img.fits", (20.0, "arcmin"), 1024)

    file_answer_testing("EVENTS", "test_events.fits", answer_store, 
                        answer_dir)
    file_answer_testing("SPECTRUM", "test_spec.fits", answer_store,
                        answer_dir)

    n_photons2, n_cells2 = make_photons("photons2", sphere, redshift, A,
                                        exp_time, thermal_model)
    n_events2 = project_photons("photons2", "events2", [1.0, -0.5, 0.2],
                                [30., 45.], absorb_model="tbabs", nH=0.1,
                                prng=prng)

    merge_files(["photons1.h5", "photons2.h5"], "photons.h5")
    merge_files(["events1.h5", "events2.h5"], "events.h5")

    with h5py.File("photons.h5", "r") as f, h5py.File("photons1.h5", "r") as f1, \
        h5py.File("photons2.h5", "r") as f2:
        assert f["data"]["energy"].size == n_photons1 + n_photons2
        assert f["data"]["x"].size == n_cells1 + n_cells2
        for k in f["data"]:
            assert_array_equal(f["data"][k][()],
                               np.concatenate([f1["data"][k][()],
                                               f2["data"][k][()]]))

    with h5py.File("events.h5", "r") as f, h5py.File("events1.h5", "r") as f1, \
            h5py.File("events2.h5", "r") as f2:
        assert f["data"]["eobs"].size == n_events1 + n_events2
        for k in f["data"]:
            assert_array_equal(f["data"][k][()],
                               np.concatenate([f1["data"][k][()],
                                               f2["data"][k][()]]))