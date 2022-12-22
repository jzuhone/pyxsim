from glob import glob

from numpy.random import RandomState
from pytest import importorskip

from pyxsim import CIESourceModel
from pyxsim.tests.utils import hdf5_answer_testing


@importorskip("yt_astro_analysis")
def test_light_cone(answer_store, answer_dir):

    from pyxsim.light_cone import XrayLightCone

    prng = RandomState(0x4D3D3D3)

    A = 20000.0
    exp_time = 1.0e6
    fov = (0.5, "deg")

    lc = XrayLightCone("enzo_tiny_cosmology/32Mpc_32.enzo", "Enzo", 0.0, 0.1, seed=24)

    source_model = CIESourceModel("apec", 0.1, 10.0, 1000, 0.3, prng=prng)

    lc.generate_events(
        "lc_photons",
        "lc_events",
        A,
        exp_time,
        fov,
        source_model,
        (30.0, 45.0),
        absorb_model="wabs",
        nH=0.02,
        sigma_pos=0.5,
        prng=prng,
    )

    fn_ph = glob("lc_photons.*.h5")
    fn_ev = glob("lc_events.*.h5")

    for fn in fn_ph:
        hdf5_answer_testing(fn, answer_store, answer_dir)
    for fn in fn_ev:
        hdf5_answer_testing(fn, answer_store, answer_dir)
