from numpy.random import RandomState
from yt.utilities.answer_testing.framework import requires_ds, \
    GenericArrayTest, data_dir_load
from pyxsim import XrayLightCone, \
    ThermalSourceModel

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

etc = "enzo_tiny_cosmology/DD0046/DD0046"

def return_data(data):
    def _return_data(name):
        return data
    return _return_data

@requires_ds(etc)
def test_light_cone():

    prng = RandomState(0x4d3d3d3)

    ds = data_dir_load(etc)

    A = 2000.
    exp_time = 1.0e5
    fov = (0.5, "deg")

    lc = XrayLightCone('%s/32Mpc_32.enzo' % etc[:-14], 'Enzo', 0., 0.1,
                       seed=24)

    source_model = ThermalSourceModel("apec", 0.1, 10.0, 1000, prng=prng)

    events = lc.generate_events(A, exp_time, fov, source_model, (30.0, 45.0),
                                absorb_model="wabs", nH=0.02,
                                sigma_pos=0.5, prng=prng)

    return_events = return_data(events.events)

    test = GenericArrayTest(ds, return_events, args=["events"])
    test_light_cone.__name__ = test.description
    yield test
