# B2 — Explicit WGS84 geodetic to local ENU transformation

Status: implementation candidate; manual Colab validation pending.

## Contract

GeoCore owns execution of geospatial coordinate transformations. This stage
adds an explicit WGS84 path:

1. geodetic latitude/longitude/ellipsoidal height to ECEF;
2. ECEF translation relative to a caller-supplied geodetic origin;
3. ECEF-to-ENU rotation at that origin.

Input order is latitude degrees, longitude degrees, ellipsoidal height metres.
Output order is East, North, Up in metres. The origin, datum, source axis order,
source units, target frame identifier, and versioned transform identifier are
returned with the result. Output arrays are defensive copies and read-only.

## Guardrails

- Only WGS84 is supported; unsupported or missing datum declarations fail.
- Origin is mandatory and must be one valid geodetic coordinate.
- No automatic axis-order, CRS, height-datum, or origin inference.
- Height is interpreted as ellipsoidal height in metres.
- This is a local tangent ENU frame; large-area use requires an explicit
  strategy for frame extent and distortion.
- Core spatial references remain descriptive and do not execute transforms.
- No pyproj dependency is introduced; the implementation uses NumPy.

## Validation

Manual Colab validation must verify exact branch/commit identity, compile,
run GeoCore's complete test suite, and run the new transform tests. No dataset
download, training, installation, GitHub Actions, push to main, or merge is
part of this candidate stage.
