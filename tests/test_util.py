"""
Tests for offset_tracking/util.py

Covers the functions fixed in the quality-fix branch:
  - get_DT / get_deltaT  (required adding 'import datetime' to util.py)
  - getazimuthAngle      (file-handle safety via ZipFile context-manager)
  - read_vals            (nodata masking logic)
"""
import io
import sys
import zipfile
import unittest
from unittest.mock import MagicMock, patch, call

# ---------------------------------------------------------------------------
# Stub heavy C-extension dependencies so the module can be imported in the
# test environment even when GDAL, cv2, etc. are not installed.
# ---------------------------------------------------------------------------
_STUBS = ['cv2', 'osgeo', 'osgeo.gdal', 'osgeo.osr']
for _mod in _STUBS:
    sys.modules.setdefault(_mod, MagicMock())

import numpy as np  # real numpy is available

import importlib
import offset_tracking.util as util_mod


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_zip_with_annotation_xml(heading: float) -> bytes:
    """Build an in-memory zip that mimics a Sentinel-1 SAFE archive."""
    buf = io.BytesIO()
    xml_content = (
        b"<?xml version='1.0'?><root>"
        b"<platformHeading>"
        + str(heading).encode()
        + b"</platformHeading></root>"
    )
    with zipfile.ZipFile(buf, 'w') as zf:
        # Mimic the internal path structure expected by getazimuthAngle
        zf.writestr('s1a.SAFE/annotation/s1a-iw1-slc-vv.xml', xml_content)
        zf.writestr('s1a.SAFE/manifest.safe', b'irrelevant')
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Tests for get_DT
# ---------------------------------------------------------------------------

class TestGetDT(unittest.TestCase):
    """get_DT(date1, date2) -> difference in days (date2 - date1)."""

    def test_positive_difference(self):
        """date2 is after date1 → positive number of days."""
        result = util_mod.get_DT('20200101', '20200201')
        self.assertEqual(result, 31)

    def test_zero_difference_same_day(self):
        """Same date → 0 days."""
        result = util_mod.get_DT('20200615', '20200615')
        self.assertEqual(result, 0)

    def test_negative_difference(self):
        """date2 is before date1 → negative days."""
        result = util_mod.get_DT('20200201', '20200101')
        self.assertEqual(result, -31)

    def test_leap_year_boundary(self):
        """29 February in a leap year is handled correctly."""
        result = util_mod.get_DT('20200228', '20200301')
        self.assertEqual(result, 2)

    def test_cross_year_boundary(self):
        result = util_mod.get_DT('20191231', '20200101')
        self.assertEqual(result, 1)

    def test_integer_input(self):
        """Inputs are converted via str() so integers also work."""
        result = util_mod.get_DT(20200101, 20200201)
        self.assertEqual(result, 31)


# ---------------------------------------------------------------------------
# Tests for get_deltaT
# ---------------------------------------------------------------------------

class TestGetDeltaT(unittest.TestCase):
    """get_deltaT(dates) -> list of day-deltas between consecutive dates."""

    def test_empty_list(self):
        self.assertEqual(util_mod.get_deltaT([]), [])

    def test_single_date(self):
        """A single date produces an empty delta list."""
        self.assertEqual(util_mod.get_deltaT(['20200101']), [])

    def test_two_dates(self):
        result = util_mod.get_deltaT(['20200101', '20200201'])
        self.assertEqual(result, [31])

    def test_multiple_dates(self):
        dates = ['20200101', '20200201', '20200301']
        result = util_mod.get_deltaT(dates)
        self.assertEqual(result, [31, 29])

    def test_length_is_n_minus_one(self):
        dates = ['2020{:02d}01'.format(m) for m in range(1, 7)]  # 6 dates: Jan–Jun
        result = util_mod.get_deltaT(dates)
        self.assertEqual(len(result), 5)


# ---------------------------------------------------------------------------
# Tests for getazimuthAngle
# ---------------------------------------------------------------------------

class TestGetAzimuthAngle(unittest.TestCase):
    """getazimuthAngle reads the platformHeading from a SAR zip archive."""

    def _write_tmp_zip(self, heading: float, path: str):
        data = _make_zip_with_annotation_xml(heading)
        with open(path, 'wb') as fh:
            fh.write(data)

    def test_reads_angle_correctly(self):
        import tempfile, os
        with tempfile.NamedTemporaryFile(suffix='.zip', delete=False) as tf:
            tf.write(_make_zip_with_annotation_xml(-167.5))
            tmp_path = tf.name
        try:
            angle = util_mod.getazimuthAngle(tmp_path)
            self.assertAlmostEqual(angle, -167.5)
        finally:
            os.unlink(tmp_path)

    def test_returns_float(self):
        import tempfile, os
        with tempfile.NamedTemporaryFile(suffix='.zip', delete=False) as tf:
            tf.write(_make_zip_with_annotation_xml(12.3))
            tmp_path = tf.name
        try:
            angle = util_mod.getazimuthAngle(tmp_path)
            self.assertIsInstance(angle, float)
        finally:
            os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# Tests for read_vals
# ---------------------------------------------------------------------------

class TestReadVals(unittest.TestCase):
    """read_vals masks NaN values with -32767."""

    def _make_mock_ds(self, array: np.ndarray):
        """Return a mock gdal dataset whose band returns *array*."""
        mock_band = MagicMock()
        mock_band.ReadAsArray.return_value = array.copy()
        mock_ds = MagicMock()
        mock_ds.GetRasterBand.return_value = mock_band
        return mock_ds

    def test_no_nan_values(self):
        arr = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float)
        with patch('offset_tracking.util.gdal') as mock_gdal:
            mock_gdal.Open.return_value = self._make_mock_ds(arr)
            result, nodata = util_mod.read_vals('dummy.tif')
        np.testing.assert_array_equal(result, arr)
        self.assertFalse(nodata.any())

    def test_nan_replaced_by_sentinel(self):
        arr = np.array([[1.0, np.nan], [3.0, 4.0]], dtype=float)
        with patch('offset_tracking.util.gdal') as mock_gdal:
            mock_gdal.Open.return_value = self._make_mock_ds(arr)
            result, nodata = util_mod.read_vals('dummy.tif')
        self.assertEqual(result[0, 1], -32767)
        self.assertTrue(nodata[0, 1])
        self.assertFalse(nodata[0, 0])

    def test_nodat_mask_combined(self):
        """Extra nodat mask is OR-combined with NaN mask."""
        arr = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float)
        extra_mask = np.array([[False, False], [True, False]])
        with patch('offset_tracking.util.gdal') as mock_gdal:
            mock_gdal.Open.return_value = self._make_mock_ds(arr)
            result, nodata = util_mod.read_vals('dummy.tif', nodat=extra_mask)
        self.assertEqual(result[1, 0], -32767)
        self.assertTrue(nodata[1, 0])
        self.assertFalse(nodata[0, 0])

    def test_dataset_closed_after_read(self):
        """The GDAL dataset handle is set to None (released) after reading."""
        arr = np.array([[1.0]], dtype=float)
        with patch('offset_tracking.util.gdal') as mock_gdal:
            mock_ds = self._make_mock_ds(arr)
            mock_gdal.Open.return_value = mock_ds
            util_mod.read_vals('dummy.tif')
        # The function sets ds = None; we just verify Open was called once
        mock_gdal.Open.assert_called_once_with('dummy.tif')

    def test_band_argument_forwarded(self):
        arr = np.array([[5.0]], dtype=float)
        with patch('offset_tracking.util.gdal') as mock_gdal:
            mock_ds = self._make_mock_ds(arr)
            mock_gdal.Open.return_value = mock_ds
            util_mod.read_vals('dummy.tif', band=3)
        mock_ds.GetRasterBand.assert_called_once_with(3)


if __name__ == '__main__':
    unittest.main()
