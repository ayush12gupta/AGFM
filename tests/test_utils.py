"""
Tests for utils.py

Covers:
  - horn_gradient        (pure numpy)
  - directional_slope    (pure numpy)
  - get_espg             (pure arithmetic)
  - Re-export verification: read_vals, get_DT, get_deltaT, getazimuthAngle
    are the same objects as in offset_tracking.util  (fix: removed duplicates)
"""
import sys
import unittest
from unittest.mock import MagicMock

# ---------------------------------------------------------------------------
# Stub C-extension dependencies before importing either module
# ---------------------------------------------------------------------------
for _mod in ['cv2', 'osgeo', 'osgeo.gdal', 'osgeo.osr',
             'psutil', 'pyproj', 'fiona', 'rasterio', 'rasterio.mask',
             'geopandas']:
    sys.modules.setdefault(_mod, MagicMock())

import numpy as np  # real numpy

# Make rasterio.mask.mask importable
import rasterio as _rasterio_mock
_rasterio_mock.mask = MagicMock()
sys.modules['rasterio.mask'] = _rasterio_mock.mask

import offset_tracking.util as util_mod
import utils as utils_mod


# ---------------------------------------------------------------------------
# Re-export identity checks
# ---------------------------------------------------------------------------

class TestUtilsReexports(unittest.TestCase):
    """
    After removing the duplicated functions from utils.py, the names in
    utils must refer to exactly the same callable objects as in
    offset_tracking.util.
    """

    def test_read_vals_is_same_object(self):
        self.assertIs(utils_mod.read_vals, util_mod.read_vals)

    def test_get_DT_is_same_object(self):
        self.assertIs(utils_mod.get_DT, util_mod.get_DT)

    def test_get_deltaT_is_same_object(self):
        self.assertIs(utils_mod.get_deltaT, util_mod.get_deltaT)

    def test_getazimuthAngle_is_same_object(self):
        self.assertIs(utils_mod.getazimuthAngle, util_mod.getazimuthAngle)


# ---------------------------------------------------------------------------
# get_espg
# ---------------------------------------------------------------------------

class TestGetEspg(unittest.TestCase):
    """get_espg computes UTM EPSG codes from latitude and zone."""

    def test_northern_hemisphere(self):
        # zone 43 N  → 32600 + 43 = 32643
        self.assertEqual(utils_mod.get_espg(30, 43), 32643)

    def test_southern_hemisphere(self):
        # zone 43 S  → 32600 + 43 + 100 = 32743
        self.assertEqual(utils_mod.get_espg(-10, 43), 32743)

    def test_equator_is_northern(self):
        # lat == 0 is treated as northern (condition is lat < 0)
        code = utils_mod.get_espg(0, 36)
        self.assertEqual(code, 32636)

    def test_zone_1(self):
        self.assertEqual(utils_mod.get_espg(45, 1), 32601)

    def test_zone_60(self):
        self.assertEqual(utils_mod.get_espg(60, 60), 32660)


# ---------------------------------------------------------------------------
# directional_slope
# ---------------------------------------------------------------------------

class TestDirectionalSlope(unittest.TestCase):
    """directional_slope combines slopex and slopey for a given azimuth angle."""

    def test_angle_zero_returns_slopey(self):
        # angle=0 → sin(0)=0, cos(0)=1 → result = slopey
        sx = np.array([1.0, 2.0])
        sy = np.array([3.0, 4.0])
        result = utils_mod.directional_slope(sx, sy, 0)
        np.testing.assert_allclose(result, sy)

    def test_angle_90_returns_slopex(self):
        # angle=90 → sin(90)=1, cos(90)=0 → result = slopex
        sx = np.array([1.0, 2.0])
        sy = np.array([3.0, 4.0])
        result = utils_mod.directional_slope(sx, sy, 90)
        np.testing.assert_allclose(result, sx, atol=1e-14)

    def test_zero_slopes_always_zero(self):
        sx = np.zeros(5)
        sy = np.zeros(5)
        for angle in [0, 45, 90, 135, 180]:
            result = utils_mod.directional_slope(sx, sy, angle)
            np.testing.assert_array_equal(result, np.zeros(5))


# ---------------------------------------------------------------------------
# horn_gradient
# ---------------------------------------------------------------------------

class TestHornGradient(unittest.TestCase):
    """horn_gradient(z, geo) computes Horn (1981) finite-difference gradients."""

    # geo is (xmin, xres, 0, ymax, 0, yres) – use unit spacing for simplicity
    UNIT_GEO = (0, 1.0, 0, 0, 0, 1.0)

    def test_flat_surface_zero_gradient(self):
        z = np.ones((5, 5), dtype=float)
        dz_dy, dz_dx = utils_mod.horn_gradient(z, self.UNIT_GEO)
        np.testing.assert_allclose(dz_dy[1:-1, 1:-1], 0, atol=1e-14)
        np.testing.assert_allclose(dz_dx[1:-1, 1:-1], 0, atol=1e-14)

    def test_linear_ramp_x_direction(self):
        """z = x (columns) → dx/dx interior ≈ 1/geo[5], dy/dx interior ≈ 0."""
        cols = np.arange(7, dtype=float)
        z = np.tile(cols, (7, 1))           # rows all identical
        dz_dy, dz_dx = utils_mod.horn_gradient(z, self.UNIT_GEO)
        # Interior gradient in x should be 1/geo[5] = 1
        np.testing.assert_allclose(dz_dx[1:-1, 1:-1], 1.0, atol=1e-12)
        # Interior gradient in y should be 0
        np.testing.assert_allclose(dz_dy[1:-1, 1:-1], 0.0, atol=1e-12)

    def test_linear_ramp_y_direction(self):
        """z = y (rows) → dy/dy interior ≈ 1/geo[1], dx/dy interior ≈ 0."""
        rows = np.arange(7, dtype=float).reshape(-1, 1)
        z = np.tile(rows, (1, 7))
        dz_dy, dz_dx = utils_mod.horn_gradient(z, self.UNIT_GEO)
        np.testing.assert_allclose(dz_dy[1:-1, 1:-1], 1.0, atol=1e-12)
        np.testing.assert_allclose(dz_dx[1:-1, 1:-1], 0.0, atol=1e-12)

    def test_output_shape_matches_input(self):
        z = np.random.rand(8, 9)
        dz_dy, dz_dx = utils_mod.horn_gradient(z, self.UNIT_GEO)
        self.assertEqual(dz_dy.shape, z.shape)
        self.assertEqual(dz_dx.shape, z.shape)

    def test_border_pixels_are_zero(self):
        """The border pixels of the output are always zero (not computed)."""
        z = np.random.rand(6, 6)
        dz_dy, dz_dx = utils_mod.horn_gradient(z, self.UNIT_GEO)
        # top/bottom rows
        np.testing.assert_array_equal(dz_dy[0, :], 0)
        np.testing.assert_array_equal(dz_dy[-1, :], 0)
        np.testing.assert_array_equal(dz_dx[0, :], 0)
        np.testing.assert_array_equal(dz_dx[-1, :], 0)
        # left/right cols
        np.testing.assert_array_equal(dz_dy[:, 0], 0)
        np.testing.assert_array_equal(dz_dy[:, -1], 0)
        np.testing.assert_array_equal(dz_dx[:, 0], 0)
        np.testing.assert_array_equal(dz_dx[:, -1], 0)


if __name__ == '__main__':
    unittest.main()
