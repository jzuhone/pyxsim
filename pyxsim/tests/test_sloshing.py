"""
Answer test pyxsim.
"""

from pyxsim import \
    ThermalSourceModel, \
    PhotonList, merge_files, EventList
from yt.utilities.answer_testing.framework import requires_ds, \
    GenericArrayTest, data_dir_load
from numpy.testing import assert_array_equal, \
    assert_allclose
from numpy.random import RandomState
from yt.units.yt_array import uconcatenate
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
    photons1 = PhotonList.from_data_source(sphere, redshift, A, exp_time,
                                           thermal_model)

    return_photons = return_data(photons1.photons)

    nphots = 0
    for i in range(4):
        phots = PhotonList.from_data_source(sphere, redshift, A, 0.25*exp_time,
                                            thermal_model)

        phots.write_h5_file("split_photons_%d.h5" % i)
        nphots += len(phots.photons["energy"])

    merge_files(["split_photons_%d.h5" % i for i in range(4)],
                "merged_photons.h5", add_exposure_times=True,
                overwrite=True)

    merged_photons = PhotonList.from_file("merged_photons.h5")
    assert len(merged_photons.photons["energy"]) == nphots
    assert merged_photons.parameters["fid_exp_time"] == exp_time

    events1 = photons1.project_photons([1.0,-0.5,0.2], [30., 45.],
                                       absorb_model="tbabs", nH=0.1, prng=prng)

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

    photons1.write_h5_file("test_photons.h5")
    events1.write_h5_file("test_events.h5")
    events1.write_fits_file("test_events.fits", 20.0, 1024)

    photons2 = PhotonList.from_file("test_photons.h5")
    events2 = EventList.from_h5_file("test_events.h5")
    events3 = EventList.from_fits_file("test_events.fits")

    for k in photons1.keys():
        if k == "energy":
            arr1 = uconcatenate(photons1[k])
            arr2 = uconcatenate(photons2[k])
        else:
            arr1 = photons1[k]
            arr2 = photons2[k]
        assert_array_equal(arr1, arr2)
    for k in events2.keys():
        assert_array_equal(events1[k], events2[k])
        assert_allclose(events2[k], events3[k], rtol=1.0e-6)

    nevents = 0

    for i in range(4):
        events = photons1.project_photons([1.0, -0.5, 0.2], [30., 45.],
                                          absorb_model="tbabs", nH=0.1,
                                          prng=prng)
        events.write_h5_file("split_events_%d.h5" % i)
        nevents += len(events["xsky"])

    merge_files(["split_events_%d.h5" % i for i in range(4)],
                "merged_events.h5", add_exposure_times=True,
                overwrite=True)

    merged_events = EventList.from_h5_file("merged_events.h5")
    assert len(merged_events["xsky"]) == nevents
    assert merged_events.parameters["exp_time"] == 4.0*exp_time

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
