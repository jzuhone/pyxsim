"""
Answer test pyxsim.
"""

from pyxsim import \
    ThermalSourceModel, \
    EventList, make_photons, \
    project_photons
from pyxsim.tests.utils import hdf5_answer_testing, file_answer_testing
from numpy.random import RandomState
import os
import tempfile
import shutil
import yt

gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"



def test_sloshing(answer_store, answer_dir):

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = RandomState(0x4d3d3d3)

    ds = yt.load(gslr)
    A = 2000.
    exp_time = 1.0e4
    redshift = 0.1

    sphere = ds.sphere("c", (0.5, "Mpc"))
    sphere.set_field_parameter("X_H", 0.75)

    thermal_model = ThermalSourceModel("apec", 0.1, 11.0, 10000, Zmet=0.3,
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

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
