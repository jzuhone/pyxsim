"""
Answer test pyxsim.
"""

from pyxsim import \
    ThermalSourceModel, \
    EventList, make_photons, \
    project_photons
from yt.utilities.answer_testing.framework import requires_ds, \
    GenericArrayTest, data_dir_load
from numpy.testing import assert_array_equal, \
    assert_allclose
from numpy.random import RandomState
import os
import tempfile
import shutil
import astropy.io.fits as pyfits


def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"


def return_data(data):
    def _return_data(name):
        return data
    return _return_data


@requires_ds(gslr)
def test_sloshing():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = RandomState(0x4d3d3d3)

    ds = data_dir_load(gslr)
    A = 2000.
    exp_time = 1.0e4
    redshift = 0.1

    sphere = ds.sphere("c", (0.1, "Mpc"))
    sphere.set_field_parameter("X_H", 0.75)

    thermal_model = ThermalSourceModel("apec", 0.1, 11.0, 10000, Zmet=0.3,
                                       thermal_broad=False, prng=prng)
    n_photons1, n_cells1 = make_photons("photons1", sphere, redshift, A, 
                                        exp_time, thermal_model)

    return_photons = return_data(photons1.photons)

    n_events1 = project_photons("photons1", "events1", [1.0, -0.5, 0.2],
                                [30., 45.], absorb_model="tbabs", nH=0.1,
                                prng=prng)

    return_events = return_data(events1.events)

    events1.write_spectrum("test_events_spec.fits", 0.2, 10.0, 2000)

    f = pyfits.open("test_events_spec.fits")
    return_spec = return_data(f["SPECTRUM"].data["COUNTS"])
    f.close()

    events1.write_fits_image("test_events_img.fits", (20.0, "arcmin"), 
                             1024)

    f = pyfits.open("test_events_img.fits")
    return_img = return_data(f[0].data)
    f.close()

    tests = [GenericArrayTest(ds, return_photons, args=["photons"]),
             GenericArrayTest(ds, return_events, args=["events"]),
             GenericArrayTest(ds, return_spec, args=["spec"]),
             GenericArrayTest(ds, return_img, args=["img"])]

    for test in tests:
        test_sloshing.__name__ = test.description
        yield test

    events1.write_fits_file("test_events.fits", 20.0, 1024)

    photons2 = PhotonList.from_file("test_photons.h5")
    events2 = EventList.from_h5_file("test_events.h5")
    events3 = EventList.from_fits_file("test_events.fits")

    for k in photons1.keys():
        arr1 = photons1[k]
        arr2 = photons2[k]
        assert_array_equal(arr1, arr2)
    for k in events2.keys():
        assert_array_equal(events1[k], events2[k])
        assert_allclose(events2[k], events3[k], rtol=1.0e-6)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
