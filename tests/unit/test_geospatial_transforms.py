import numpy as np
import pytest

from vyra_geocore.geospatial.transforms import (
    GeospatialTransformError,
    ecef_to_enu,
    geodetic_to_ecef,
    geodetic_to_enu,
)


def test_origin_transforms_to_zero_enu():
    origin = (37.0, -122.0, 12.0)
    result = geodetic_to_enu(
        np.array([origin]), origin_geodetic=origin, datum="WGS84"
    )
    np.testing.assert_allclose(result.coordinates, np.zeros((1, 3)), atol=1e-8)
    assert result.target_axis_order == "east_north_up"
    assert result.target_units == "metres"
    assert result.transform_id.endswith("_v1")
    assert not result.coordinates.flags.writeable


def test_ecef_known_equator_point():
    result = geodetic_to_ecef(np.array([0.0, 0.0, 0.0]))
    np.testing.assert_allclose(result, [6378137.0, 0.0, 0.0], atol=1e-6)


def test_east_displacement_has_positive_east_component():
    origin = (0.0, 0.0, 0.0)
    point = np.array([[0.0, 0.001, 0.0]])
    enu = geodetic_to_enu(point, origin_geodetic=origin, datum="WGS84").coordinates
    assert enu[0, 0] > 100.0
    assert abs(enu[0, 1]) < 1.0


def test_rejects_implicit_or_unsupported_datum_and_axis_order():
    with pytest.raises(GeospatialTransformError):
        geodetic_to_enu(np.zeros((1, 3)), origin_geodetic=(0, 0, 0), datum="")
    with pytest.raises(GeospatialTransformError):
        geodetic_to_enu(
            np.zeros((1, 3)), origin_geodetic=(0, 0, 0),
            datum="WGS84", axis_order="longitude_latitude_height"
        )


@pytest.mark.parametrize("bad", [
    [91.0, 0.0, 0.0],
    [0.0, 181.0, 0.0],
    [0.0, 0.0, float("nan")],
])
def test_rejects_invalid_geodetic_coordinates(bad):
    with pytest.raises(GeospatialTransformError):
        geodetic_to_ecef(np.array(bad))


def test_ecef_to_enu_requires_explicit_single_origin():
    with pytest.raises(GeospatialTransformError):
        ecef_to_enu(np.zeros((1, 3)), np.zeros((2, 3)))
