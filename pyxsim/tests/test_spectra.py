import soxs
from numpy.testing import assert_allclose

from pyxsim.spectral_models import IGMSpectralModel, TableCIEModel


def test_apec():

    amod = TableCIEModel("apec", 0.1, 10.0, 10000, 1.0, 10.0, thermal_broad=True)
    amod.prepare_spectrum(0.2)

    acspec, amspec, _ = amod.get_spectrum(6.0)
    spec = acspec + 0.3 * amspec

    agen = soxs.ApecGenerator(0.1, 10.0, 10000, broadening=True)
    spec2 = agen.get_spectrum(6.0, 0.3, 0.2, 1.0e-14)

    assert_allclose(spec[0, :], spec2.flux.value * spec2.de.value)

    emin = 0.5
    emax = 7.0

    pf = amod.make_fluxf(emin, emax, energy=False)
    ef = amod.make_fluxf(emin, emax, energy=True)

    eidxs = (amod.ebins[:-1] > emin) & (amod.ebins[1:] < emax)
    c, m, _ = pf(6.0)
    ec, em, _ = ef(6.0)

    assert_allclose(c + 0.3 * m, (spec2.flux * spec2.de)[eidxs].value.sum())
    assert_allclose(
        ec + 0.3 * em, (spec2.emid * spec2.flux * spec2.de)[eidxs].value.sum()
    )


def test_igm():

    imod = IGMSpectralModel(0.2, 3.0, 1000)
    imod.prepare_spectrum(0.05)

    icspec, imspec, _ = imod.get_spectrum(1.0, 0.01)
    spec = icspec + 0.3 * imspec

    igen = soxs.IGMGenerator(0.2, 3.0, 1000)
    spec2 = igen.get_spectrum(1.0, 0.01, 0.3, 0.05, 1.0e-14)

    assert_allclose(spec[0, :], spec2.flux.value * spec2.de.value)

    icspec, imspec, _ = imod.get_spectrum(2.0, 0.01)
    spec3 = icspec + 0.3 * imspec

    cgen = soxs.CloudyCIEGenerator(0.2, 3.0, 1000)
    spec4 = cgen.get_spectrum(2.0, 0.3, 0.05, 1.0e-14)

    assert_allclose(spec3[0, :], spec4.flux.value * spec4.de.value)

    emin = 0.5
    emax = 2.0

    pf = imod.make_fluxf(emin, emax, energy=False)
    ef = imod.make_fluxf(emin, emax, energy=True)

    eidxs = (imod.ebins[:-1] > emin) & (imod.ebins[1:] < emax)
    c, m, _ = pf(1.0, 0.01)
    ec, em, _ = ef(1.0, 0.01)

    assert_allclose(c + 0.3 * m, (spec2.flux * spec2.de)[eidxs].value.sum())
    assert_allclose(
        ec + 0.3 * em, (spec2.emid * spec2.flux * spec2.de)[eidxs].value.sum()
    )

    c2, m2, _ = pf(2.0, 0.01)
    ec2, em2, _ = ef(2.0, 0.01)

    assert_allclose(c2 + 0.3 * m2, (spec4.flux * spec4.de)[eidxs].value.sum())
    assert_allclose(
        ec2 + 0.3 * em2, (spec4.emid * spec4.flux * spec4.de)[eidxs].value.sum()
    )
