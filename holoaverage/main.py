# This file is part of holoaverage.
# Copyright (c) 2018 Tore Niermann
#
# holoaverage is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with holoaverage.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

import sys
import os.path
import numpy as np
import json
import io

from holoaverage.defocus import propagate
from .series import DataSet, LazyLoadingSeries
from .dm3 import load_dm3
from .hdf5 import saveHDF5, loadHDF5
from .reconstruction import series_reconstruction, holo_reconstruction
from .average import holoAverage
from .rawalign import rawAlign, extractROI
from .camera import ParameterizedMTF
from .jsoncleaner import remove_trailing_commas, remove_comments
from .filter import FilterFunction

def print_syntax(prog_name):
    print('Syntax:')
    print('\t%s [-v] parameter-file' % prog_name)
    print()
    print('The parameter file is a JSON file. See documentation for details.')
    print('If "-" is passed as parameter file name, the parameters are read from stdin.')
    print()
    print('Options:')
    print('\t-v Verbose')


def parse_parameters(param):
    # Remove comments and trailing commas
    param = remove_comments(param)
    param = remove_trailing_commas(param)

    # Load as JSON
    return json.loads(param)


def load_file(filename):
    """Return loader function depending on extension"""
    parts = filename.split("?", 1)
    ext = os.path.splitext(parts[0])[1].lower()
    if ext == ".dm3":
        return load_dm3(filename)
    elif ext == ".hdf5":
        return loadHDF5(filename)
    else:
        raise ValueError("Unrecognized image extension.")


def rescale_fourier(data, shape=None, out=None):
    """
    Rescales `data` to `shape` using zero-padding/cropping in Fourier space.

    :param data: Input array
    :param shape: Output shape
    :param out: Output array
    """
    if len(data.shape) != 2:
        raise ValueError("Input data must be 2D")
    if shape is None:
        shape = out.shape
    if len(shape) != 2:
        raise ValueError("Destination shape must be 2D")
    if data.dtype.kind != 'c':
        raise ValueError("Input data expected to be complex.")
    half = (min(shape[0], data.shape[0]) - 1) // 2, (min(shape[1], data.shape[1]) - 1) // 2

    if out is None:
        out = np.empty(shape, dtype=data.dtype)
    elif out.shape != shape:
        raise ValueError("Output array shape must fit given shape.")
    out.fill(0)
    out[shape[0] // 2 - half[0]:shape[0] // 2 + half[0], shape[1] // 2 - half[1]:shape[1] // 2 + half[1]] = \
        np.fft.fftshift(np.fft.fft2(data))[data.shape[0] // 2 - half[0]:data.shape[0] // 2 + half[0],
                                           data.shape[1] // 2 - half[1]:data.shape[1] // 2 + half[1]]
    out[...] = np.fft.ifft2(np.fft.ifftshift(out))
    return out


def create_synthetic_empty(camera_dev, sideband_pos, empty_size, sampling, mask_type):
    # Synthesize empty
    if camera_dev[0].shape != camera_dev[1].shape:
        raise ValueError("X and Y distortion arrays must have same shape")

    empty = DataSet(camera_dev[0].shape, dtype=np.float32)
    empty.attrs["dim_scale"] = sampling
    empty.attrs["dim_unit"] = ["nm"] * 2
    carrier_px = float(sideband_pos[0] - empty.shape[1] // 2) / empty.shape[1], \
                 float(sideband_pos[1] - empty.shape[0] // 2) / empty.shape[0]
    rx = np.arange(0, empty.shape[1]).reshape(1, -1)
    ry = np.arange(0, empty.shape[0]).reshape(-1, 1)
    empty.array[...] = np.cos(2.0 * np.pi * ((rx + camera_dev[0]) * carrier_px[0] + (ry + camera_dev[1]) * carrier_px[1]))

    carrier = carrier_px[0] / sampling[0], carrier_px[1] / sampling[1]
    empty_holo = holo_reconstruction(empty, carrier=carrier, shape=(empty_size, empty_size), maskType=mask_type)
    return empty_holo


def reconstruct_average(param, basepath="", verbose=0):
    """
    Reconstruct averaged holograms

    :param param: Dictionary with parameters
    :param basepath: All filenames are taken relative to this path.
    :param verbose: Verbosity level
    """
    # Get required parameters
    if 'object_names' in param:
        object_names = str(param['object_names'])
        object_first = int(param.get('object_first'))
        object_last = int(param.get('object_last'))
    else:
        object_names = None
    if 'object_size' in param:
        object_size = int(param['object_size'])
    else:
        object_size = None

    synthesize_empty = bool(param.get('synthesize_empty', False))
    if synthesize_empty:
        empty_names = None
        empty_override = None
    elif 'empty_override' in param:
        empty_override = str(param['empty_override'])
        empty_names = None
    elif 'empty_names' in param:
        empty_names = str(param['empty_names'])
        empty_first = int(param['empty_first'])
        empty_last = int(param['empty_last'])
        empty_override = None
    else:
        empty_names = None
        empty_override = None
    if 'empty_size' in param:
        empty_size = int(param['empty_size'])
    else:
        empty_size = None
    sideband_pos = np.array(param['sideband_pos'], dtype=float)
    if 'roi' in param:
        roi = np.array(param['roi'], dtype=int)
    else:
        roi = None
    cut_off = float(param['cut_off'])

    # Get optional parameters
    path = os.path.abspath(os.path.join(basepath, param.get('path', '')))
    if 'sampling' in param:
        sampling = np.ones(2, dtype=float) * param['sampling']
    else:
        sampling = None
    if 'voltage' in param:
        voltage = float(param["voltage"])
    else:
        voltage = None
    if 'binning' in param:
        binning = int(param["binning"])
    else:
        binning = None
    filter_func = param.get('filter_func', 'edge')
    only_phase = bool(param.get('only_phase', False))
    if 'align_roi' in param:
        align_roi = param.get('align_roi')
        if align_roi is not None:
            align_roi = np.array(align_roi, dtype=int)
            disable_raw_align = False
        else:
            disable_raw_align = True
    else:
        align_roi = None
        disable_raw_align = False
    adjust_defocus = bool(param.get('adjust_defocus', False))
    adjust_shift = bool(param.get('adjust_shift', True))
    adjust_tilt = bool(param.get('adjust_tilt', False))
    output_series = bool(param.get('output_series', False))
    output_aligned = bool(param.get('output_aligned', False))
    output_file = param.get('output')
    if 'output' in param:
        output_file = os.path.join(path, str(param['output']))
    defocus_first = float(param.get('defocus_first', 0.0))
    defocus_step = float(param.get('defocus_step', 0.0))
    empty_exclude = np.array(param.get('empty_exclude', []), dtype=int)
    object_exclude = np.array(param.get('object_exclude', []), dtype=int)
    if 'mtf' in param:
        mtf = ParameterizedMTF(param["mtf"])
    else:
        mtf = None
    distortions = param.get('camera_distortions')

    # Adjust filter
    mask_type = FilterFunction(max_q=cut_off, mask_type=filter_func)

    # Reconstruct empty wave
    if synthesize_empty:
        # Create synthetic empty
        camera_dev_names = [os.path.join(path, name) for name in distortions]
        camera_dev = [load_file(name).array for name in camera_dev_names]
        if verbose > 0:
            print("Synthesizing empty hologram from\n\t%s\n\t%s" % (camera_dev_names[0], camera_dev_names[1]))
        if empty_size is None:
            if object_size is None:
                raise ValueError("Either parameter 'empty_size' or parameter 'object_size' are needed.")
            empty_size = object_size
        if sampling is None:
            # Get sampling from first object file
            if object_names is None:
                raise ValueError("Either parameter 'sampling' or parameter 'object_names' is needed.")
            object_index = [index for index in range(object_first, object_last + 1) if index not in object_exclude]
            object_files = [os.path.join(path, object_names % index) for index in object_index]
            empty_sampling = np.empty(2, dtype=float)
            empty_sampling[...] = load_file(object_files[0]).attrs["dim_scale"]
        else:
            empty_sampling = np.ones(2, dtype=float) * sampling

        empty = create_synthetic_empty(camera_dev, sideband_pos, empty_size, empty_sampling, mask_type)
        if output_file:
            saveHDF5(output_file, empty, dataName="empty")
    elif empty_override is not None:
        empty_override = os.path.join(path, empty_override)
        empty = load_file(empty_override)
    elif empty_names is not None:
        empty_files = [os.path.join(path, empty_names % index) for index in range(empty_first, empty_last + 1) if index not in empty_exclude]
        empty_series = LazyLoadingSeries.fromFiles(empty_files, load_file, verbose=1)
        if sampling is not None:
            empty_series.attrs['dim_scale'] = sampling
            empty_series.attrs['dim_unit'] = ['nm'] * 2
            empty_series.grid = None
        if voltage is not None:
            empty_series.attrs["voltage(kV)"] = voltage
        if binning is not None:
            empty_series.attrs["binning"] = (binning, binning)
        carrier = float(sideband_pos[0] - empty_series.shape[1] // 2) / empty_series.shape[1] / empty_series.attrs['dim_scale'][0], \
                  float(sideband_pos[1] - empty_series.shape[0] // 2) / empty_series.shape[0] / empty_series.attrs['dim_scale'][1]

        if empty_size is None:
            raise ValueError("Parameter 'empty_size' is missing.")
        empty_holo = series_reconstruction(empty_series, carrier=carrier, shape=(empty_size, empty_size), verbose=verbose, mtf=mtf, maskType=mask_type)
        del empty_series

        if len(empty_holo) > 1:
            empty = holoAverage(empty_holo, adjustDefocus=False, adjustShift=False, adjustTilt=adjust_tilt, iterations=4, verbose=verbose * 2)
        else:
            empty = empty_holo[0]

        if output_file:
            saveHDF5(output_file, empty, dataName="empty")
    else:
        if empty_size is None:
            if object_size is None:
                raise ValueError("Either parameter 'empty_size' or parameter 'object_size' are needed.")
            empty_size = object_size
        empty = DataSet(shape=(empty_size, empty_size), dtype=np.complex64)
        empty.array.fill(1.0)

    if object_names is None:
        return

    # Align object wave
    object_index = [index for index in range(object_first, object_last + 1) if index not in object_exclude]
    object_files = [os.path.join(path, object_names % index) for index in object_index]
    data_series = LazyLoadingSeries.fromFiles(object_files, load_file, verbose=1)
    if sampling is not None:
        data_series.attrs['dim_scale'] = np.ones(2, dtype=float) * sampling
        data_series.attrs['dim_unit'] = ['nm'] * 2
        data_series.grid = None
    if voltage is not None:
        data_series.attrs["voltage(kV)"] = voltage
    if binning is not None:
        data_series.attrs["binning"] = (binning, binning)
    data_shape = data_series.shape
    carrier = float(sideband_pos[0] - data_shape[1] // 2) / data_shape[1] / data_series.attrs['dim_scale'][0], \
              float(sideband_pos[1] - data_shape[0] // 2) / data_shape[0] / data_series.attrs['dim_scale'][1]

    if roi is None:
        roi = [0, 0, data_series.shape[1], data_series.shape[0]]
    data_series.attrs['roi'] = roi
    if align_roi is None:
        align_roi = roi
    if not disable_raw_align and len(data_series) > 1:
        data_series = rawAlign(data_series, qMax=cut_off, roi=align_roi, verbose=verbose)
        data_series.attrs['align_roi'] = align_roi
    data_rois = extractROI(data_series, roi, verbose=verbose)
    del data_series

    # Save series
    if output_aligned and output_file:
        data_rois.saveHDF5(output_file, 'aligned_rois')

    # Reconstruct object series
    if object_size is None:
        raise ValueError("Parameter 'object_size' is missing.")
    holo_series = series_reconstruction(data_rois, carrier=carrier, shape=(object_size, object_size), verbose=verbose, mtf=mtf, maskType=mask_type)
    del data_rois

    # Zeropad empty
    big_empty = rescale_fourier(empty.array, data_shape)
    if only_phase:
        big_empty = np.exp(1j * np.angle(big_empty))

    # Correct holo_series with empty hologram
    roi_shape = roi[3] - roi[1], roi[2] - roi[0]
    big_data = np.zeros(roi_shape, dtype=np.complex64)
    for flat_index in range(holo_series.indexSize):
        index = np.unravel_index(flat_index, holo_series.indexShape)

        # Zeropad data
        rescale_fourier(holo_series[index].array, out=big_data)

        tmp_roi = holo_series[index].attrs.get("roi", roi)
        sy, sx = tmp_roi[1], tmp_roi[0]
        ny, nx = roi_shape
        dy, dx = 0, 0
        if sy < 0:
            ny += sy
            dy -= sy
            sy = 0
        if sx < 0:
            nx += sx
            dx -= sx
            sx = 0
        if (sy + ny) > data_shape[0]:
            ny = data_shape[0] - sy
        if (sx + nx) > data_shape[1]:
            nx = data_shape[1] - sx

        big_data[dy:dy + ny, dx:dx + nx] /= big_empty[sy:sy + ny, sx:sx + nx]
        rescale_fourier(big_data, out=holo_series[index].array)

    # Save series
    if output_series and output_file:
        holo_series.saveHDF5(output_file, 'series')

    # Average series
    defocus = np.array([(index - object_first) * defocus_step + defocus_first for index in object_index])
    if len(holo_series) > 1:
        out, var = holoAverage(holo_series, defocus=defocus, adjustDefocus=adjust_defocus, adjustShift=adjust_shift, adjustTilt=adjust_tilt, verbose=verbose * 2, variance=True)
    else:
        out = propagate(holo_series[0], defocus[0])
        var = None
    if output_file:
        saveHDF5(output_file, out)
    if output_file and (var is not None):
        saveHDF5(output_file, var, dataName="variance")


def main(argv=sys.argv):
    """
    Main program
    :param argv: Arguments
    :returns: Exit code
    """
    arg_index = 1
    param_file = None
    verbose = 0
    while arg_index < len(argv):
        if argv[arg_index][0] == '-':
            option = argv[arg_index][1:]
            arg_index += 1
            if option == 'v':
                verbose += 1
            elif option == '':
                param_file = '-'
            else:
                print("Unknown option:", option)
                print()
                print_syntax()
                return 2    # Use exit code 2 for Syntax
        else:
            param_file = argv[arg_index]
            arg_index += 1

    # Missing parameters ?
    if param_file is None:
        print("Parameter file missing.")
        print()
        print_syntax(argv[0])
        return 2  # Use exit code 2 for Syntax

    # Read param
    if param_file == '-':
        if verbose > 0:
            print("Loading parameters from stdin")
        param = sys.stdin.read()
        basepath = os.path.abspath(os.getcwd())
    else:
        param_file = os.path.abspath(param_file)
        if verbose > 0:
            print("Loading parameters from\n\t%s" % param_file)
        with io.open(param_file, "rt", encoding="utf-8") as file:
            param = file.read()
        basepath = os.path.dirname(param_file)
    param = parse_parameters(param)

    # Reconstruct
    reconstruct_average(param, basepath, verbose=verbose)
    return 0  # Use exit code 0 for success


