from pyxsim import TableCIEModel
from numpy.testing import assert_allclose


def test_apec():

    amod = TableCIEModel("apec", 0.1, 10.0, 10000, thermal_broad=True)
    amod.prepare_spectrum(0.2, 1, 10)

    acspec, amspec, _ = amod.get_spectrum(6.0)
    spec = acspec+0.3*amspec

    spec2 = amod.return_spectrum(6.0, 0.3, 0.2, 1.0e-14)

    assert_allclose(spec[0,:], spec2.v)
