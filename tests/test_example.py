import unittest
import os
import tempfile
import h5py
import numpy as np

from holoaverage.main import reconstruct_average, rescale_fourier

# Setup path to example data
EXAMPLE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../example/GaN_holographic_focal_series'))

# MTF for example data
EXAMPLE_MTF = [["CONSTANT", -2.25536738e-02],
               ["LORENTZIAN", 1.02543658e-05, 1.15367655e-04],
               ["LORENTZIAN", 2.49224357e-02, 5.35262063e-02],
               ["GAUSSIAN", 4.60461599e-01, 4.36842560e+02]]

class TextExamples(unittest.TestCase):
    """Test whether the examples run."""

    GAN_DATA_CONVERGENCE_256 = 1.337721e8           # Expected error for 256px object reconstruction
    GAN_DATA_CONVERGENCE_256_FOCUS = 1.249292e8     # Expected error for 256px object reconstruction with defocus adjustment
    GAN_EMPTY_CONVERGENCE_384 = 5.681450e10         # Expected error for 384px empty reconstruction

    # Test parameters
    VERBOSE = 0
    DELETE_OUTPUT = True

    def setUp(self):
        self.temp_outputs = []

    def tearDown(self):
        if self.DELETE_OUTPUT:
            while len(self.temp_outputs) > 0:
                filename = self.temp_outputs.pop()
                os.unlink(filename)

    def touch_temp_output(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".hdf5") as tmpfile:
             filename = tmpfile.name
        self.temp_outputs.append(filename)
        return filename

    def default_GaN_param(self):
        # Touch temporary file
        # Setup parameter structure
        param = {
            "object_names": os.path.join(EXAMPLE_PATH, "a1_%03d.dm3"),
            "object_first": 1,
            "object_last": 20,
            "object_size": 256,
            "empty_names": os.path.join(EXAMPLE_PATH, "empty_%03d.dm3"),
            "empty_first": 1,
            "empty_last": 20,
            "empty_size": 384,
            "sampling": 0.00519824,
            "defocus_first": 20.0,
            "defocus_step": -2.0,
            "sideband_pos" : [1136, 1304],
            "roi" : [128, 128, 1920, 1920],
            "output" : self.touch_temp_output(),
            "filter_func" : ["BUTTERWORTH", 14],
            "cut_off" : 12.0,
            "mtf": EXAMPLE_MTF
        }
        return param

    def test_GaN_example(self):
        param = self.default_GaN_param()
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            data_convergence = output["data"].attrs["convergence"]
            empty_convergence = output["empty"].attrs["convergence"]

        # Test for empty convergence
        self.assertTrue(np.allclose(empty_convergence[-1], self.GAN_EMPTY_CONVERGENCE_384, rtol=1e-4))
        delta = (empty_convergence[1:] - empty_convergence[:-1]) / empty_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

        # Test for data convergence
        self.assertTrue(np.allclose(data_convergence[-1], self.GAN_DATA_CONVERGENCE_256, rtol=1e-4))
        delta = (data_convergence[1:] - data_convergence[:-1]) / data_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

    def test_GaN_example_other_sizes(self):
        param = self.default_GaN_param()
        param["object_size"] = 384
        param["empty_size"] = 512
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            data_convergence = output["data"].attrs["convergence"]
            empty_convergence = output["empty"].attrs["convergence"]

        # Test for empty convergence
        self.assertTrue(np.allclose(empty_convergence[-1], self.GAN_EMPTY_CONVERGENCE_384 * (384.0/512.0)**2, rtol=0.1))
        delta = (empty_convergence[1:] - empty_convergence[:-1]) / empty_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

        # Test for data convergence
        self.assertTrue(np.allclose(data_convergence[-1], self.GAN_DATA_CONVERGENCE_256 * (256.0/384.0)**2, rtol=0.1))
        delta = (data_convergence[1:] - data_convergence[:-1]) / data_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

    def test_GaN_example_disabled_pyfftw(self):
        from holoaverage.fft import pyfftw_present, disable_pyfftw
        if not pyfftw_present():
            return      # Disable test: already tested by other test
        disable_pyfftw()

        param = self.default_GaN_param()
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            data_convergence = output["data"].attrs["convergence"]
            empty_convergence = output["empty"].attrs["convergence"]

        # Test for empty convergence
        self.assertTrue(np.allclose(empty_convergence[-1], self.GAN_EMPTY_CONVERGENCE_384, rtol=1e-4))
        delta = (empty_convergence[1:] - empty_convergence[:-1]) / empty_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

        # Test for data convergence
        self.assertTrue(np.allclose(data_convergence[-1], 1.337719e8, rtol=1e-4))
        delta = (data_convergence[1:] - data_convergence[:-1]) / data_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

    def test_GaN_focus_and_tilt_alignment(self):
        param = self.default_GaN_param()
        param["adjust_defocus"] = True
        param["adjust_tilt"] = True
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            data_convergence = output["data"].attrs["convergence"]
            empty_convergence = output["empty"].attrs["convergence"]

        # Test for empty convergence
        self.assertTrue(np.allclose(empty_convergence[-1], self.GAN_EMPTY_CONVERGENCE_384, rtol=1e-1))
        delta = (empty_convergence[1:] - empty_convergence[:-1]) / empty_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

        # Test for data convergence
        self.assertTrue(np.allclose(data_convergence[-1], self.GAN_DATA_CONVERGENCE_256_FOCUS, rtol=1e-4))
        delta = (data_convergence[1:] - data_convergence[:-1]) / data_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

    def test_GaN_no_mtf(self):
        param = self.default_GaN_param()
        param["object_last"] = 5    # Smaller series for speed
        param["empty_last"] = 5
        del param["mtf"]
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            data_convergence = output["data"].attrs["convergence"]
            empty_convergence = output["empty"].attrs["convergence"]

        # Test for empty convergence
        delta = (empty_convergence[1:] - empty_convergence[:-1]) / empty_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

        # Test for data convergence
        delta = (data_convergence[1:] - data_convergence[:-1]) / data_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

    def test_GaN_no_objects(self):
        param = self.default_GaN_param()
        del param["object_names"]
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            self.assertFalse("data" in output)
            self.assertFalse("variance" in output)
            empty_convergence = output["empty"].attrs["convergence"]

        # Test for empty convergence
        self.assertTrue(np.allclose(empty_convergence[-1], self.GAN_EMPTY_CONVERGENCE_384, rtol=1e-4))
        delta = (empty_convergence[1:] - empty_convergence[:-1]) / empty_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

    def test_GaN_no_empty(self):
        param = self.default_GaN_param()
        del param["empty_names"]
        del param["empty_size"]
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            self.assertFalse("empty" in output)
            data_convergence = output["data"].attrs["convergence"]

        # Test for data convergence
        self.assertTrue(np.allclose(data_convergence[-1], 5.6592e14, rtol=1e-4))
        delta = (data_convergence[1:] - data_convergence[:-1]) / data_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

    def test_GaN_separate_empty_and_object(self):
        # Reconstruct empty
        empty_param = self.default_GaN_param()
        del empty_param["object_names"]
        reconstruct_average(empty_param, verbose=self.VERBOSE)
        with h5py.File(empty_param["output"], "r") as output:
            empty_convergence = output["empty"].attrs["convergence"]

        # Test for empty convergence
        self.assertTrue(np.allclose(empty_convergence[-1], self.GAN_EMPTY_CONVERGENCE_384, rtol=1e-4))
        delta = (empty_convergence[1:] - empty_convergence[:-1]) / empty_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

        # Reconstruct data
        object_param = self.default_GaN_param()
        object_param["empty_names"] = "/dev/null/some_invalid_filename"
        object_param["empty_override"] = empty_param["output"] + "?empty"
        reconstruct_average(object_param, verbose=self.VERBOSE)
        with h5py.File(object_param["output"], "r") as output:
            data_convergence = output["data"].attrs["convergence"]

        # Test for data convergence
        self.assertTrue(np.allclose(data_convergence[-1], self.GAN_DATA_CONVERGENCE_256, rtol=1e-4))
        delta = (data_convergence[1:] - data_convergence[:-1]) / data_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

    def test_GaN_single_empty_and_object(self):
        # Reconstruct empty
        param = self.default_GaN_param()
        param["object_last"] = 1
        param["empty_last"] = 1
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            data = output["data"][...]
            empty = output["empty"][...]

        # Test for shape
        self.assertEqual(data.shape, (param["object_size"],) * 2)
        self.assertEqual(empty.shape, (param["empty_size"],) * 2)

    def test_GaN_full_size(self):
        # Reconstruct empty
        param = self.default_GaN_param()
        param["object_last"] = 1
        param["empty_last"] = 1
        param["object_size"] = 2048
        param["empty_size"] = 2048
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            data = output["data"][...]
            empty = output["empty"][...]

        # Test for shape
        self.assertEqual(data.shape, (param["object_size"],) * 2)
        self.assertEqual(empty.shape, (param["empty_size"],) * 2)

    def test_GaN_synthetic_empty(self):
        # Reconstruct empty for displacements
        empty_param = self.default_GaN_param()
        del empty_param["object_names"]
        carrier = (np.array(empty_param["sideband_pos"], dtype=float) - 1024) / 2048
        carrier /= np.dot(carrier, carrier)
        reconstruct_average(empty_param, verbose=self.VERBOSE)
        with h5py.File(empty_param["output"], "r") as output:
            empty = output["empty"][...]

        # Create displacements from empty hologram
        shift = np.angle(rescale_fourier(empty, (2048, 2048))) / 2.0 / np.pi
        displacement_file = self.touch_temp_output()
        with h5py.File(displacement_file, "w") as output:
            output.create_dataset("dx", data=shift * carrier[0])
            output.create_dataset("dy", data=shift * carrier[1])

        # Reconstruct data
        object_param = self.default_GaN_param()
        object_param["empty_names"] = "/dev/null/some_invalid_filename"
        object_param["empty_override"] = "/dev/null/some_invalid_filename"
        object_param["synthesize_empty"] = True
        object_param["camera_distortions"] = [displacement_file + "?dx", displacement_file + "?dy"]
        reconstruct_average(object_param, verbose=self.VERBOSE)
        with h5py.File(object_param["output"], "r") as output:
            new_empty = output["empty"][...]
            data_convergence = output["data"].attrs["convergence"]

        # Test for empty equivalence
        np.testing.assert_allclose(np.angle(empty / new_empty)[10:-10, 10:-10], 0.0, atol=0.1)

        # Test for data convergence
        self.assertTrue(np.allclose(data_convergence[-1], 5.1516e+11, rtol=1e-4))
        delta = (data_convergence[1:] - data_convergence[:-1]) / data_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

    def test_GaN_non_square_roi(self):
        param = self.default_GaN_param()
        param["object_last"] = 5    # Smaller series for speed
        param["empty_last"] = 5
        param["roi"] = [256, 640, 1792, 1408]
        reconstruct_average(param, verbose=self.VERBOSE)
        with h5py.File(param["output"], "r") as output:
            data_convergence = output["data"].attrs["convergence"]
            empty_convergence = output["empty"].attrs["convergence"]

        # Test for empty convergence
        delta = (empty_convergence[1:] - empty_convergence[:-1]) / empty_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))

        # Test for data convergence
        delta = (data_convergence[1:] - data_convergence[:-1]) / data_convergence[:-1]
        self.assertTrue(np.all(delta < +1e-4))


if __name__ == '__main__':
    unittest.main()

