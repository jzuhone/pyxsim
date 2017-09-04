from pyxsim import TableApecModel
from yt.utilities.answer_testing.framework import \
    GenericArrayTest
from yt.testing import requires_module, fake_random_ds
from numpy.testing import assert_allclose

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

ds = fake_random_ds(64)
@requires_module("astropy")
def test_apec():

    amod = TableApecModel(0.1, 10.0, 10000, thermal_broad=True)
    amod.prepare_spectrum(0.2)

    acspec, amspec, _ = amod.get_spectrum(6.0)
    spec = acspec+0.3*amspec

    spec2 = amod.return_spectrum(6.0, 0.3, 0.2, 1.0e-14)

    assert_allclose(spec.v, spec2.v)

    def spec_test():
        return spec.v

    test = GenericArrayTest(ds, spec_test)
    test_apec.__name__ = test.description
    yield test
