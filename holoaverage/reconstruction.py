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
import warnings

from .series import Series, DataSet
from .filter import FilterFunction
from .grid import Grid
from .fft import empty_aligned, zeros_aligned, fromHermite


__all__ = ("HoloReconstructor", "holo_reconstruction", "series_reconstruction")


class HoloReconstructor(object):
    """
    Class that handles reconstruction of holograms.

    The parameters passed to the __init__ function are used to reconstruct the
    individual holograms, passed to the apply function.
    """

    # Members:
    #    _grid        Grid object
    #    _mask        Mask
    #    _destScale   Scale of destination data
    #    _destShape   Scale of destination data
    #    _destRect    Indices (dx, dy, dx+nx, dy+ny) of destination region in FFT-space
    #    _srcRect     Indices (sx, sy, sx+nx, sy+ny) of source region in FFT-space
    #    _carrier     Carrier frequencies (x,y) in 1/nm
    #    _destIFFT    IFFT for destination data
    def __init__(self, grid, mask=None, carrier=None, shape=None, mtf=None):
        """
        Arguments:
            grid
                Grid object for the data to be reconstructed.
                If a DataSet is passed, Grid.fromDataSet is used.
            shape
                Size of reconstruction (YDim, XDim) (Defaults to 1/4 size of series shape)
            mask
                If mask is an instance of FilterFunction, that instance is used.
                If mask is a reciprocal frequency in 1/nm, i.e. a number,  this frequency as cutoff is used.
                Otherwise a mask with 1/3 carrier frequency is used.
                In the latter two cases the type of the mask is the default.
            carrier
                Carrier frequency in [1/nm]
            mtf
                mtf function (optional: same shape as grid)
        """
        # Parse parameters
        if not isinstance(grid, Grid):
            self._grid = Grid.fromDataSet(grid)
        else:
            self._grid = grid
        if len(self._grid.shape) != 2:
            raise ValueError("2D series expected.")
        if shape is None:
            shape = tuple(s // 4 for s in self._grid.shape)
        # Calculate scales
        src_scale = grid.realSampling
        # Test for pixels out of desired region
        carrier_px = np.dot(src_scale, np.array((carrier[1], carrier[0]), dtype=float)) * np.array(self._grid.shape, dtype=float)
        ry, rx = self._grid.shape[0] // 2, self._grid.shape[1] // 2
        ny, nx = shape
        hy, hx = ny // 2, nx // 2
        dy, dx = 0, 0
        sy = ry + int(carrier_px[0]) - hy
        sx = rx + int(carrier_px[1]) - hx
        # print sx, sy, rx, ry, carrier_px, hx, hy
        if sy < 0:
            ny += sy
            dy -= sy
            sy = 0
        if sx < 0:
            nx += sx
            dx -= sx
            sx = 0
        if (sy + ny) > self._grid.shape[0]:
            ny = self._grid.shape[0] - sy
        if (sx + nx) > self._grid.shape[1]:
            nx = self._grid.shape[1] - sx
        # Masks
        if mask is None:
            mask = FilterFunction(max_q=np.sqrt(np.sum(carrier ** 2)) / 3.0)
        elif not isinstance(mask, FilterFunction):
            mask = FilterFunction(max_q=float(mask))
        self._mask = mask.calculate_shifted(self._grid, roi=[sx, sy, sx + nx, sy + ny], origin=carrier)
        # MTF
        if mtf is not None:
            self._mask /= np.fft.fftshift(mtf)[sy:sy + ny, sx:sx + nx]
        # Save parameters
        self._destScale = src_scale * np.array(self._grid.shape) / np.array(shape)
        self._destShape = shape
        self._carrier = carrier
        self._destRect = (dx, dy, dx + nx, dy + ny)
        self._srcRect = (sx, sy, sx + nx, sy + ny)
        self._reco_attrs = {'dim_scale': self._destScale,
                            'carrier(nm-1)': self._carrier,
                            'space': +1,
                            'reconstructionCutOff2(nm-2)': mask.max_q,
                            'reconstructionMaskType': mask.mask_type}
        # Prepare FFT
        self._reco_grid = Grid(shape, realSampling=self._destScale, dtype=self._grid.complexType)
        self._reco_grid.prepareFFT(forwardFFT=False, backwardFFT=True, forwardRFFT=False, backwardRFFT=False)

    def apply(self, data):
        if (data.shape != self._grid.shape):
            raise ValueError("'data' has wrong shape.")
        sub = empty_aligned(self._grid.shape, dtype=self._grid.floatType)
        side = zeros_aligned(self._destShape, dtype=self._grid.complexType)
        sub[...] = data.array
        if data.attrs.get("space", 0) >= 0:
            tmp = np.fft.fftshift(fromHermite(self._grid.forwardRFFT(sub), sub.shape))
        else:
            tmp = sub
        # Extract and iFFT
        dx0, dy0, dx1, dy1 = self._destRect
        sx0, sy0, sx1, sy1 = self._srcRect
        side[dy0:dy1, dx0:dx1] = tmp[sy0:sy1, sx0:sx1] * self._mask
        result = DataSet(self._destShape, dtype=self._grid.complexType)
        result.array[...] = self._reco_grid.backwardFFT(np.fft.ifftshift(side))
        result.attrs.update(data.attrs)
        result.attrs.update(self._reco_attrs)
        return result

    def __call__(self, data):
        return self.apply(data)

    @property
    def grid(self):
        return self._grid

    @property
    def shape(self):
        return self._destShape

    @property
    def carrier(self):
        return self._carrier


def holo_reconstruction(data, carrier=None, qMax=None, qMax2=None, shape=None, mtf=None, maskType=None):
    """
    Do a holographic reconstruction of the data and returns the result as new dataset.

    Arguments:
        data
            The dataset to be reconstructed
        shape
            Size of reconstruction (YDim, XDim) (Defaults to 1/4 size of series shape)
        carrier
            Carrier frequency in [1/nm]
        qMax
            highest spatial frequency considered in [1/nm] (Defaults to 1/3 carrier frequency)
        qMax2
            highest spatial frequency squared considered in [1/nm^2] (Defaults to 1/3 carrier frequency)
            overrides qMax
        mtf
            MTF (see :class:`ParameterizedMTF`)
        maskType
            "edge", tuple("butterworth", order), "gaussian" (uses qMax as 1/e**2)
            Alternatively an instance of FilterFunction
    """
    grid = Grid.fromDataSet(data)
    # MTF
    if mtf is not None:
        mtf = mtf.mtf(grid.shape, data.attrs.get("binning", 1), dtype=grid.floatType)
    else:
        mtf = None
    if not isinstance(maskType, FilterFunction):
        if maskType is not None:
            warnings.warn("Using an explicit mask type is deprecated.", DeprecationWarning)
        if qMax2 is not None:
            qMax = None
        maskType = FilterFunction(max_q=qMax, max_q2=qMax2, mask_type=maskType)
    reconstructor = HoloReconstructor(grid, carrier=carrier, shape=shape, mask=maskType, mtf=mtf)
    return reconstructor(data)


def series_reconstruction(series, carrier=None, qMax=None, qMax2=None, shape=None, mtf=None, verbose=0, maskType=None):
    """Do a holographic reconstruction of the data and returns
    the result as a new series.

    Arguments:
        series
            The series to be reconstructed
        shape
            Size of reconstruction (YDim, XDim) (Defaults to 1/4 size of series shape)
        carrier
            Carrier frequency in [1/nm]
        qMax
            highest spatial frequency considered in [1/nm] (Defaults to 1/3 carrier frequency)
        qMax2
            highest spatial frequency squared considered in [1/nm^2] (Defaults to 1/3 carrier frequency)
            overrides qMax
        mtf
            MTF (see :class:`ParameterizedMTF`)
        maskType
            "EDGE", tuple("BUTTERWORTH", order), "GAUSSIAN" (uses qMax as 1/e**2)
            Alternatively an instance of FilterFunction
    """
    grid = series.grid
    # MTF
    if mtf is not None:
        mtf = mtf.mtf(grid.shape, series.attrs.get("binning", 1), dtype=grid.floatType)
    else:
        mtf = None
    if not isinstance(maskType, FilterFunction):
        warnings.warn("Using an explicit mask type is deprecated.", DeprecationWarning)
        if qMax2 is not None:
            qMax = None
        maskType = FilterFunction(max_q=qMax, max_q2=qMax2, mask_type=maskType)
    reconstructor = HoloReconstructor(grid, carrier=carrier, shape=shape, mask=maskType, mtf=mtf)
    destSeries = Series(series.indexShape, reconstructor.shape, dtype=grid.complexType)
    if verbose > 0:
        print("Reconstructing...")
        print("\t", end="")
    for flatIndex in range(series.indexSize):
        index = np.unravel_index(flatIndex, series.indexShape)
        if verbose > 0:
            print(". ", end="")
        destSeries[index] = reconstructor(series[index])
    if verbose > 0:
        print()
    index0 = np.unravel_index(0, series.indexShape)
    destSeries.attrs.update(series.attrs)
    destSeries.attrs['dim_scale'] = destSeries[index0].attrs['dim_scale']
    destSeries.attrs['carrier(nm-1)'] = destSeries[index0].attrs['carrier(nm-1)']
    destSeries.attrs['space'] = destSeries[index0].attrs['space']
    destSeries.attrs['reconstructionCutOff2(nm-2)'] = destSeries[index0].attrs['reconstructionCutOff2(nm-2)']
    destSeries.attrs['reconstructionMaskType'] = destSeries[index0].attrs['reconstructionMaskType']
    return destSeries
