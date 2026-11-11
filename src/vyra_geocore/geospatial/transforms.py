"""Explicit WGS84 geodetic -> ECEF -> local ENU transformations.

Coordinates use [latitude_deg, longitude_deg, ellipsoidal_height_m].
ENU output uses [east_m, north_m, up_m] relative to an explicit origin.
No CRS guessing, implicit origin, or axis-order inference is performed.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

WGS84_A_M = 6378137.0
WGS84_INV_F = 298.257223563
WGS84_F = 1.0 / WGS84_INV_F
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
TRANSFORM_ID = "wgs84_geodetic_to_ecef_to_local_enu_v1"


class GeospatialTransformError(ValueError):
    """Raised when geodetic coordinates or transform metadata are invalid."""


@dataclass(frozen=True)
class ENUTransformResult:
    coordinates: np.ndarray
    source_datum: str
    source_axis_order: str
    source_units: str
    target_frame: str
    target_axis_order: str
    target_units: str
    origin_geodetic: tuple[float, float, float]
    transform_id: str = TRANSFORM_ID

    def __post_init__(self) -> None:
        values = np.asarray(self.coordinates, dtype=np.float64)
        if values.ndim < 1 or values.shape[-1] != 3 or not np.all(np.isfinite(values)):
            raise GeospatialTransformError("ENU coordinates must be finite (..., 3) values.")
        values = values.copy()
        values.setflags(write=False)
        object.__setattr__(self, "coordinates", values)


def _validate_geodetic(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim < 1 or arr.shape[-1] != 3:
        raise GeospatialTransformError(f"{name} must have shape (..., 3).")
    if not np.all(np.isfinite(arr)):
        raise GeospatialTransformError(f"{name} contains non-finite values.")
    lat, lon = arr[..., 0], arr[..., 1]
    if np.any((lat < -90.0) | (lat > 90.0)):
        raise GeospatialTransformError("Latitude must be within [-90, 90] degrees.")
    if np.any((lon < -180.0) | (lon > 180.0)):
        raise GeospatialTransformError("Longitude must be within [-180, 180] degrees.")
    return arr


def geodetic_to_ecef(coordinates: np.ndarray) -> np.ndarray:
    """Convert WGS84 [lat_deg, lon_deg, ellipsoidal_height_m] to ECEF metres."""
    geo = _validate_geodetic(coordinates, "coordinates")
    lat = np.deg2rad(geo[..., 0])
    lon = np.deg2rad(geo[..., 1])
    height = geo[..., 2]
    sin_lat, cos_lat = np.sin(lat), np.cos(lat)
    sin_lon, cos_lon = np.sin(lon), np.cos(lon)
    prime_vertical = WGS84_A_M / np.sqrt(1.0 - WGS84_E2 * sin_lat**2)
    x = (prime_vertical + height) * cos_lat * cos_lon
    y = (prime_vertical + height) * cos_lat * sin_lon
    z = (prime_vertical * (1.0 - WGS84_E2) + height) * sin_lat
    return np.stack((x, y, z), axis=-1)


def ecef_to_enu(ecef_coordinates: np.ndarray, origin_geodetic: np.ndarray) -> np.ndarray:
    """Rotate ECEF metre coordinates into ENU at an explicit WGS84 origin."""
    ecef = np.asarray(ecef_coordinates, dtype=np.float64)
    if ecef.ndim < 1 or ecef.shape[-1] != 3 or not np.all(np.isfinite(ecef)):
        raise GeospatialTransformError("ECEF coordinates must be finite (..., 3) values.")
    origin = _validate_geodetic(np.asarray(origin_geodetic, dtype=np.float64), "origin")
    if origin.shape != (3,):
        raise GeospatialTransformError("origin_geodetic must be exactly one [lat, lon, h] tuple.")
    origin_ecef = geodetic_to_ecef(origin)
    delta = ecef - origin_ecef
    lat0, lon0 = np.deg2rad(origin[:2])
    slat, clat = np.sin(lat0), np.cos(lat0)
    slon, clon = np.sin(lon0), np.cos(lon0)
    rotation = np.array([
        [-slon, clon, 0.0],
        [-slat * clon, -slat * slon, clat],
        [clat * clon, clat * slon, slat],
    ], dtype=np.float64)
    return np.einsum("ij,...j->...i", rotation, delta)


def geodetic_to_enu(
    coordinates: np.ndarray,
    *,
    origin_geodetic: tuple[float, float, float],
    datum: str,
    axis_order: str = "latitude_longitude_height",
    units: str = "degrees_degrees_metres",
    target_frame: str = "local_enu",
) -> ENUTransformResult:
    """Perform the complete, explicitly declared WGS84 geodetic -> ECEF -> ENU transform."""
    if not isinstance(datum, str) or not datum.strip():
        raise GeospatialTransformError("datum must be an explicit non-empty string.")
    if datum.strip().upper() != "WGS84":
        raise GeospatialTransformError("Only explicitly declared WGS84 datum is supported.")
    if axis_order != "latitude_longitude_height":
        raise GeospatialTransformError("axis_order must be 'latitude_longitude_height'.")
    if units != "degrees_degrees_metres":
        raise GeospatialTransformError("units must be 'degrees_degrees_metres'.")
    if not isinstance(target_frame, str) or not target_frame.strip():
        raise GeospatialTransformError("target_frame must be a non-empty explicit identifier.")
    origin = _validate_geodetic(np.asarray(origin_geodetic, dtype=np.float64), "origin")
    if origin.shape != (3,):
        raise GeospatialTransformError("origin_geodetic must be one [lat, lon, h] tuple.")
    ecef = geodetic_to_ecef(coordinates)
    enu = ecef_to_enu(ecef, origin)
    return ENUTransformResult(
        coordinates=enu,
        source_datum="WGS84",
        source_axis_order=axis_order,
        source_units=units,
        target_frame=target_frame.strip(),
        target_axis_order="east_north_up",
        target_units="metres",
        origin_geodetic=tuple(float(v) for v in origin),
    )
