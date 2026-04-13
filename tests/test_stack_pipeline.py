"""
Tests for stack_pipeline.py  →  file_pair_selection()

Covers the fix: file handles are now opened with 'with' statements so they
are always properly closed even when an exception occurs.

Also covers the pairing logic: ascending / descending filenames whose
acquisition dates are within *num* days of each other are matched.
"""
import io
import sys
import unittest
from unittest.mock import MagicMock, patch, mock_open, call

# ---------------------------------------------------------------------------
# Stub heavy dependencies so stack_pipeline can be imported
# ---------------------------------------------------------------------------
for _mod in ['stack_process', 'batch3Dinversion', 'pandas']:
    sys.modules.setdefault(_mod, MagicMock())

# Make numpy available as a real module (used below in assertions)
import numpy as _np  # real numpy is installed

from stack_pipeline import file_pair_selection


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_filename(date_str: str, orbit: str = 'IW') -> str:
    """
    Build a minimal fake Sentinel-1 filename where token [5] (0-indexed)
    starts with an 8-digit date.

    Pattern: A_B_C_D_E_YYYYMMDDTHHMMSS_...
    """
    return f'S1A_IW_SLC__1SDV_{date_str}T000000_IW_{orbit}.zip'


def _write_txt(filenames):
    """Produce newline-separated file content (with trailing newline)."""
    return '\n'.join(filenames) + '\n'


class FilePairSelectionTest(unittest.TestCase):

    # ------------------------------------------------------------------
    # File-handle safety
    # ------------------------------------------------------------------

    def test_file_handles_closed_after_normal_run(self):
        """
        Both file handles must be closed after a normal (non-error) call.
        We verify this by using mock_open and checking that __exit__ was
        called (which is what 'with' guarantees).
        """
        asc_content = _write_txt([_make_filename('20200101')])
        des_content = _write_txt([_make_filename('20200101')])

        asc_mock = mock_open(read_data=asc_content)
        des_mock = mock_open(read_data=des_content)

        with patch('builtins.open', side_effect=[asc_mock.return_value,
                                                  des_mock.return_value]):
            file_pair_selection('asc.txt', 'des.txt')

        # __exit__ is called when the 'with' block exits normally
        asc_mock.return_value.__exit__.assert_called_once()
        des_mock.return_value.__exit__.assert_called_once()

    def test_file_handles_closed_on_exception(self):
        """
        Even if processing raises an exception, the context manager
        guarantees __exit__ is called (file is closed).
        We simulate this by injecting a bad filename that causes a parse
        error after opening, and confirm __exit__ was still invoked.
        """
        # An ascending file with a malformed name (can't split token [5])
        bad_content = 'BADFILENAME\n'
        good_content = _write_txt([_make_filename('20200101')])

        asc_mock = mock_open(read_data=bad_content)
        des_mock = mock_open(read_data=good_content)

        with patch('builtins.open', side_effect=[asc_mock.return_value,
                                                  des_mock.return_value]):
            with self.assertRaises(Exception):
                file_pair_selection('asc.txt', 'des.txt')

        # The 'with open(...)' block must have been entered → __exit__ called
        asc_mock.return_value.__exit__.assert_called_once()

    # ------------------------------------------------------------------
    # Pairing logic
    # ------------------------------------------------------------------

    def _run(self, asc_names, des_names, num=4):
        asc_content = _write_txt(asc_names)
        des_content = _write_txt(des_names)
        with patch('builtins.open',
                   side_effect=[mock_open(read_data=asc_content).return_value,
                                 mock_open(read_data=des_content).return_value]):
            return file_pair_selection('asc.txt', 'des.txt', num=num)

    def test_exact_same_date_is_paired(self):
        """
        Regression test: the original code had no branch for diff == 0,
        causing an infinite loop when ascending and descending dates matched
        exactly.  The fix uses `if diff <= 0` so exact matches are always
        treated as within-range and both indices are advanced.
        """
        asc = [_make_filename('20200101')]
        des = [_make_filename('20200101')]
        asc_fl, des_fl = self._run(asc, des)
        self.assertEqual(asc_fl, asc)
        self.assertEqual(des_fl, des)

    def test_within_num_days_is_paired(self):
        asc = [_make_filename('20200101')]
        des = [_make_filename('20200103')]  # 2 days apart, num=4
        asc_fl, des_fl = self._run(asc, des, num=4)
        self.assertEqual(len(asc_fl), 1)
        self.assertEqual(len(des_fl), 1)

    def test_beyond_num_days_not_paired(self):
        asc = [_make_filename('20200101')]
        des = [_make_filename('20200110')]  # 9 days apart, num=4
        asc_fl, des_fl = self._run(asc, des, num=4)
        self.assertEqual(asc_fl, [])
        self.assertEqual(des_fl, [])

    def test_exactly_at_boundary_is_paired(self):
        """A difference of exactly *num* days should be paired (abs(diff)>num is False)."""
        asc = [_make_filename('20200101')]
        des = [_make_filename('20200105')]  # diff = 4, num = 4 → paired
        asc_fl, des_fl = self._run(asc, des, num=4)
        self.assertEqual(len(asc_fl), 1)

    def test_one_beyond_boundary_not_paired(self):
        """A difference of num+1 days should NOT be paired."""
        asc = [_make_filename('20200101')]
        des = [_make_filename('20200106')]  # diff = 5, num = 4 → skipped
        asc_fl, des_fl = self._run(asc, des, num=4)
        self.assertEqual(asc_fl, [])

    def test_multiple_pairs_all_matched(self):
        asc = [_make_filename('20200101'), _make_filename('20200201')]
        des = [_make_filename('20200101'), _make_filename('20200201')]
        asc_fl, des_fl = self._run(asc, des)
        self.assertEqual(len(asc_fl), 2)
        self.assertEqual(len(des_fl), 2)

    def test_unpaired_excess_files_are_skipped(self):
        """If one track has extra files with no matching partner, they are dropped."""
        asc = [_make_filename('20200101'), _make_filename('20200201'),
               _make_filename('20200301')]
        des = [_make_filename('20200101')]          # only one partner
        asc_fl, des_fl = self._run(asc, des)
        self.assertEqual(len(asc_fl), 1)
        self.assertEqual(len(des_fl), 1)

    def test_empty_ascending_returns_empty(self):
        # Empty files produce an empty split list [''], which has one element;
        # but the date parsing would fail. Use a truly empty pair of lists.
        asc_content = ''
        des_content = ''
        # split('\n')[:-1] on '' returns [] so this is safe
        with patch('builtins.open',
                   side_effect=[mock_open(read_data=asc_content).return_value,
                                 mock_open(read_data=des_content).return_value]):
            asc_fl, des_fl = file_pair_selection('asc.txt', 'des.txt')
        self.assertEqual(asc_fl, [])
        self.assertEqual(des_fl, [])

    def test_result_lengths_match(self):
        """Ascending and descending result lists must always have the same length."""
        asc = [_make_filename(f'2020{m:02d}01') for m in range(1, 7)]
        des = [_make_filename(f'2020{m:02d}01') for m in range(1, 7)]
        asc_fl, des_fl = self._run(asc, des)
        self.assertEqual(len(asc_fl), len(des_fl))


if __name__ == '__main__':
    unittest.main()
