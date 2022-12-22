import astropy.wcs as pywcs
import numpy as np
from numpy.random import RandomState
from numpy.testing import assert_allclose
from yt import YTArray, YTQuantity

from pyxsim.lib.sky_functions import pixel_to_cel


def test_pixel_to_cel():

    prng = RandomState(24)

    n_evt = 100000

    sky_center = YTArray([30.0, 45.0], "deg")

    rr = YTQuantity(100.0, "kpc") * prng.uniform(size=n_evt)
    theta = 2.0 * np.pi * prng.uniform(size=n_evt)
    xx = rr * np.cos(theta)
    yy = rr * np.sin(theta)

    D_A = YTQuantity(100.0, "Mpc")

    d_a = D_A.to("kpc").v

    xx = xx.d / d_a
    yy = yy.d / d_a

    xsky1 = xx.copy()
    ysky1 = yy.copy()

    pixel_to_cel(xsky1, ysky1, sky_center.d)

    xx = np.rad2deg(xx) * 3600.0  # to arcsec
    yy = np.rad2deg(yy) * 3600.0  # to arcsec

    # We set a dummy pixel size of 1 arcsec just to compute a WCS
    dtheta = 1.0 / 3600.0

    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = [0.0, 0.0]
    wcs.wcs.crval = list(sky_center)
    wcs.wcs.cdelt = [-dtheta, dtheta]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = ["deg"] * 2

    xsky2, ysky2 = wcs.wcs_pix2world(xx, yy, 1)

    assert_allclose(xsky1, xsky2)
    assert_allclose(ysky1, ysky2)


if __name__ == "__main__":
    test_pixel_to_cel()
