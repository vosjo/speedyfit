import pytest
from speedyfit.photometry_query import get_parallax


def test_get_parallax():

    plx, e_plx = get_parallax('BD+34 1543', radius=3)

    assert plx == pytest.approx(5.172398, 0.000001)
    assert e_plx == pytest.approx(0.0397, 0.0001)
