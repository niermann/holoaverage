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
import numpy as np

from .grid import Grid, ScaleMatrix
from .filter import FilterFunction
from .series import Series, DataSet
from .fft import empty_aligned, zeros_aligned


def rawAlign(series, qMax=None, verbose=0, roi=None, filter=None):
    """Does a raw alignment on pixel scale of the data.
    The result is stored in the attribute raw_shift of the series.

    Arguments:
        series
            The series to be aligned
        qMax
            highest spatial frequency considered in rcpr units (alternatively use filter)
        roi
            Image region (lowerX, lowerY, ..., upperX, upperY, ...) to use. Defaults to whole dataset
        filter
            Filter function (qMax must be None)

    Returns original series with 'raw_shift' attribute
    """
    if verbose > 1:
        import matplotlib.pyplot as plt
    seriesNDim = len(series.shape)
    if verbose > 0:
        spacing = " " * ((seriesNDim - 1) * 7)
        print("Raw alignment...")
        print("\t[NN] D[px]  " + spacing + "T[px]  " + spacing)
        print(("\t[%02d]" + " %+5d" * (2 * seriesNDim)) % ((0,) + (0,) * (2 * seriesNDim)))
    # Prepare ROI
    if roi is None:
        roiShape = series.shape
        roiSlice = (slice(None),) * seriesNDim
        grid = series.grid
    else:
        roiShape = tuple((roi[seriesNDim + i] - roi[i] for i in range(seriesNDim - 1, -1, -1)))
        roiSlice = tuple((slice(roi[i], roi[seriesNDim + i]) for i in range(seriesNDim - 1, -1, -1)))
        grid = Grid(roiShape, realSampling=series.grid.realSampling, dtype=series.grid.complexType)
    if qMax is not None:
        if filter is not None:
            raise ValueError("Either 'qMax' or 'filter' must be given, not both.")
        mask = FilterFunction(max_q=qMax, mask_type="EDGE").calculate(grid)
    elif isinstance(filter, FilterFunction):
        mask = filter.calculate(grid)
    else:
        raise ValueError("Either 'qMax' or 'filter' must be given.")

    # Get shifts
    rawShift = np.zeros(series.indexShape + (seriesNDim,), dtype=int)
    sumShift = np.zeros(len(series.shape), dtype=int)
    thisFFT = empty_aligned(roiShape, dtype=series.grid.complexType)
    thisFFT[...] = series[np.unravel_index(0, series.indexShape)].array[roiSlice]
    thisFFT = grid.forwardFFT(thisFFT)
    lastFFT = empty_aligned(roiShape, dtype=series.grid.complexType)
    for i in range(1, series.indexSize):
        lastFFT, thisFFT = thisFFT, lastFFT
        thisFFT[...] = series[np.unravel_index(i, series.indexShape)].array[roiSlice]
        thisFFT = grid.forwardFFT(thisFFT)
        scale = 0.5 * abs(thisFFT[0, 0]) + 0.5 * abs(lastFFT[0, 0])
        tmp = lastFFT * np.conj(thisFFT)
        tmp /= (abs(tmp) + 0.0001 * scale)
        ccf = zeros_aligned(roiShape, dtype=series.grid.complexType)
        np.multiply(tmp, mask, ccf)
        if verbose > 1:
            plt.subplot(121)
            plt.imshow(np.fft.fftshift(np.log(abs(ccf) + 1)) + np.fft.fftshift(np.log(abs(lastFFT) + 1)), cmap="gray")
        ccf = np.fft.fftshift(grid.backwardFFT(ccf)).real
        if verbose > 1:
            plt.subplot(122)
            plt.imshow(ccf, cmap="gray", extent=[-ccf.shape[1] // 2, ccf.shape[1] // 2, -ccf.shape[0] // 2, ccf.shape[0] // 2])
        pos = np.argmax(ccf)
        pos = np.array(np.unravel_index(pos, roiShape)) - (np.array(ccf.shape) // 2)
        pos = pos[::-1]
        sumShift += pos
        rawShift[np.unravel_index(i, series.indexShape) + (Ellipsis,)] = sumShift
        if verbose > 0:
            print(("\t[%02d]" + " %+5d" * (2 * seriesNDim)) % ((i,) + tuple(pos) + tuple(sumShift)))
        if verbose > 1:
            plt.show()
    series.attrs["raw_shift"] = rawShift.copy()
    return series


def extractROI(series, roi, binning=1, verbose=0, dtype=None):
    """
    Extracts ROI from a series the result as a new series.

    Respects 'raw_shift' attribute.

    Arguments:
        series
            The series to be reconstructed
        roi
            Region of interest in [px] (left, top, right, bottom)
        binning
            Binning (X,Y) to use
        dtype
            Type of ROI series, defaults to series type
    """
    if len(series.shape) != 2:
        raise ValueError("2D series expected.")
    shape = (roi[3] - roi[1], roi[2] - roi[0])
    binning = np.diag(ScaleMatrix(2).set(binning).getMatrix()).astype(int)
    if (shape[0] % binning[1]) != 0 or (shape[1] % binning[0]) != 0:
        raise ValueError("ROI shape must be a multiple of binning.")
    if dtype is None:
        dtype = series.dtype
    raw_shift = series.attrs.get("raw_shift", np.zeros(series.indexShape + (2,), dtype=int))
    # Adjust attributes
    dim_scale = ScaleMatrix(2).set(series.attrs.get('dim_scale'))
    dim_offset = series.attrs.get('dim_offset', np.zeros(2, dtype=float))
    if dim_scale.get() is not None:
        dim_offset = np.dot(dim_scale.getMatrix(), roi[0:2]) + dim_offset
        dim_scale.set(dim_scale.getMatrix() * binning.astype(float))
    old_binning = np.diag(ScaleMatrix(2).set(series.attrs.get('binning', 1)).getMatrix())
    overrides = {"dim_scale": dim_scale.getCompact(minRank=1),
                 "dim_offset": dim_offset,
                 "binning": old_binning * binning}

    # Extracting
    result = Series(series.indexShape, (shape[0] / binning[1], shape[1] / binning[0]), dtype=dtype)
    result.attrs.update(series.attrs)
    result.attrs.update(overrides)
    result.attrs['roi'] = roi
    if verbose > 0:
        print("Extracting ROI...")
        print("\t", end="")
    for flatIndex in range(series.indexSize):
        index = np.unravel_index(flatIndex, series.indexShape)
        if verbose > 0:
            print(". ", end="")
        # Get ROI
        ny, nx = shape[0], shape[1]
        _y, _x = -raw_shift[index + (1,)], -raw_shift[index + (0,)]
        sy, sx = roi[1] + _y, roi[0] + _x
        dy, dx = 0, 0
        if sy < 0:
            ny += sy
            dy -= sy
            sy = 0
        if sx < 0:
            nx += sx
            dx -= sx
            sx = 0
        if (sy + ny) > series.shape[0]:
            ny = series.shape[0] - sy
        if (sx + nx) > series.shape[1]:
            nx = series.shape[1] - sx
        tmp = np.empty(shape, dtype=dtype)
        tmp[...] = np.mean(series[index].array)
        tmp[dy:dy + ny, dx:dx + nx] = series[index].array[sy:sy + ny, sx:sx + nx]
        data = DataSet(result.shape, dtype=dtype)
        data.array[...] = 0
        for i in range(binning[1]):
            for j in range(binning[0]):
                data.array[...] += tmp[i::binning[1], j::binning[0]]
        data.attrs.update(series[index].attrs)
        data.attrs.update(overrides)
        data.attrs['roi'] = (_x + roi[0], _y + roi[1], _x + roi[2], _y + roi[3])
        result[index] = data
    if verbose > 0:
        print()
    return result
