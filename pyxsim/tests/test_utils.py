from astropy.units import Quantity
from yt import YTQuantity

from pyxsim.utils import parse_value


def test_parse_value():
    t_yt = parse_value(YTQuantity(300.0, "ks"), "s")
    t_astropy = parse_value(Quantity(300.0, "ks"), "s")
    t_float = parse_value(300000.0, "s")
    t_tuple = parse_value((300.0, "ks"), "s")

    assert t_astropy == t_yt
    assert t_float == t_yt
    assert t_tuple == t_yt
