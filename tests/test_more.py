import unittest
import os
import tempfile
import h5py
import numpy as np

import warnings
warnings.filterwarnings("error",category=DeprecationWarning)
warnings.filterwarnings("error",category=PendingDeprecationWarning)

from holoaverage.main import holoaverage

# Setup path to example data
DATA_FILENAME = os.path.abspath(os.path.join(os.path.dirname(__file__), 'more-test-data.hdf5'))

# MTF for 200kV
MTF200 = [["CONSTANT", -1.804734e-02],
          ["LORENTZIAN", 1.762429e-02, 3.134162e-02],
          ["GAUSSIAN", 5.302573e-02, 1.327864e+04],
          ["GAUSSIAN", 3.982183e-01, 1.926439e+02]]


class TestMore(unittest.TestCase):
    """Several tests."""

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

    def touch_temp_output(self, suffix=".hdf5"):
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmpfile:
             filename = tmpfile.name
        self.temp_outputs.append(filename)
        return filename

    def default_n1_param(self):
        return {
            "voltage" : 200.0,
            "binning" : 1,
            "sampling": 1.97218,
            "object_names" : DATA_FILENAME + "?dataset=n1_0",
            "object_first" : 0,
            "object_last" : 0,
            "object_size" : 384,
            "empty_names" : DATA_FILENAME + "?dataset=n1_0",
            "empty_first" : 0,
            "empty_last" : 0,
            "empty_size" : 512,
            "sideband_pos" : [920, 1178],
            "roi" : [119, 366, 1655, 1902],
            "enable_raw_alignment" :  False,
            "filter_func": ["BUTTERWORTH", 14],
            "cut_off" : 0.008,
            "only_phase" : True,
            "adjust_shift" : False,
            "mtf" : MTF200,
            "output_name": self.touch_temp_output(),
        }

    @staticmethod
    def inner_quarter_stats(param):
        with h5py.File(param["output_name"], "r") as output:
            dataset = output[param.get("output_prefix", "") + "data"]

            qs = dataset.shape[0] // 4, dataset.shape[1] // 4
            phase_data = np.angle(dataset[qs[0]:3*qs[0], qs[1]:3*qs[1]])

            phase_mean = np.mean(phase_data)
            phase_std = np.std(phase_data)

            if TestMore.VERBOSE >= 1:
                print(phase_mean, phase_std)
            if TestMore.VERBOSE >= 2:
                import matplotlib.pyplot as plt
                plt.imshow(np.angle(dataset))
                plt.show()

        return phase_mean, phase_std

    def test_same_empty_as_object(self):
        param = self.default_n1_param()
        holoaverage(param, verbose=self.VERBOSE)
        phase_mean, phase_std = self.inner_quarter_stats(param)
        self.assertLess(abs(phase_std), 0.15)

    def test_odd_sizes(self):
        param = self.default_n1_param()

        param["object_size"] = 3**5
        holoaverage(param, verbose=self.VERBOSE)
        phase_mean, phase_std = self.inner_quarter_stats(param)
        self.assertLess(abs(phase_std), 0.15)
        param["object_size"] = 384

        param["empty_size"] = 3**5
        holoaverage(param, verbose=self.VERBOSE)
        phase_mean, phase_std = self.inner_quarter_stats(param)
        self.assertLess(abs(phase_std), 0.15)
        param["empty_size"] = 512

        param["roi"] = [param["roi"][0], param["roi"][1], param["roi"][0] + 1215, param["roi"][1] + 1215]
        holoaverage(param, verbose=self.VERBOSE)
        phase_mean, phase_std = self.inner_quarter_stats(param)
        self.assertLess(abs(phase_std), 0.15)
