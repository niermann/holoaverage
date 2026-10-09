import unittest
import os
import io
import json
import shutil
import tempfile

try:
    from unittest import mock
except ImportError:
    import mock

from holoaverage.main import holoaverage, main


class TestOutputName(unittest.TestCase):
    """Lightweight tests for the handling of the 'output_name' parameter (no reconstruction is performed)."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def write_param_file(self, param, filename="my-param.json"):
        param_file = os.path.join(self.temp_dir, filename)
        with io.open(param_file, "wt", encoding="utf-8") as file:
            file.write(json.dumps(param))
        return param_file

    def run_main(self, argv):
        """Run main() with holoaverage() mocked, returns (exit code, param passed to holoaverage)."""
        with mock.patch("holoaverage.main.holoaverage") as mocked:
            exit_code = main(["holoaverage"] + argv)
        self.assertEqual(mocked.call_count, 1)
        return exit_code, mocked.call_args[0][0]

    def test_default_from_param_file(self):
        param_file = self.write_param_file({"sideband_pos": [0, 0]})
        exit_code, param = self.run_main([param_file])
        self.assertEqual(exit_code, 0)
        self.assertEqual(param["output_name"], os.path.join(self.temp_dir, "my-param.hdf5"))

    def test_default_from_param_file_without_extension(self):
        param_file = self.write_param_file({"sideband_pos": [0, 0]}, filename="my-param")
        exit_code, param = self.run_main([param_file])
        self.assertEqual(exit_code, 0)
        self.assertEqual(param["output_name"], os.path.join(self.temp_dir, "my-param.hdf5"))

    def test_param_file_with_hdf5_extension_fails(self):
        param_file = self.write_param_file({"sideband_pos": [0, 0]}, filename="my-param.HDF5")
        with mock.patch("holoaverage.main.holoaverage") as mocked:
            exit_code = main(["holoaverage", param_file])
        self.assertEqual(exit_code, 1)
        self.assertEqual(mocked.call_count, 0)

    def test_param_file_with_hdf5_extension_and_explicit_output_name(self):
        param_file = self.write_param_file({"sideband_pos": [0, 0], "output_name": "other.hdf5"},
                                           filename="my-param.hdf5")
        exit_code, param = self.run_main([param_file])
        self.assertEqual(exit_code, 0)
        self.assertEqual(param["output_name"], "other.hdf5")

    def test_explicit_output_name_kept(self):
        param_file = self.write_param_file({"sideband_pos": [0, 0], "output_name": "other.hdf5"})
        exit_code, param = self.run_main([param_file])
        self.assertEqual(exit_code, 0)
        self.assertEqual(param["output_name"], "other.hdf5")

    def test_deprecated_output_kept(self):
        param_file = self.write_param_file({"sideband_pos": [0, 0], "output": "other.hdf5"})
        exit_code, param = self.run_main([param_file])
        self.assertEqual(exit_code, 0)
        self.assertNotIn("output_name", param)
        self.assertEqual(param["output"], "other.hdf5")

    def test_no_default_from_stdin(self):
        with mock.patch("sys.stdin", io.StringIO(u'{"sideband_pos": [0, 0]}')):
            exit_code, param = self.run_main(["-"])
        self.assertEqual(exit_code, 0)
        self.assertNotIn("output_name", param)

    def test_missing_output_name_raises(self):
        param = {"sideband_pos": [0, 0], "cut_off": 0.1}
        with self.assertRaisesRegex(ValueError, "output_name"):
            holoaverage(param)
