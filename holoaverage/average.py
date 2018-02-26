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
from math import pi, atan2
from scipy.optimize import leastsq

from .series import DataSet
from .fft import aligned_copy
from .grid import ScaleMatrix
from .defocus import calcWaveLength

__all__ = ['HoloAveraging', 'holoAverage']


class HoloAveraging(object):
    # Private members
    #    _grid           The grid object
    #    _rx, _ry        Sampling point of real grid
    #    _qx, _qy        Frequencies of reciprocal grid
    #    _qx2, _qy2      Squared frequencies of reciprocal grid
    #    _margin         Fraction of shape used as margin in comparisons
    #    _series         The original series
    #    _shift          Shifts
    #    _tilt           Tilts
    #    _defocus        Defocus
    #    _factor         Prefactor
    #    _error          Per image error
    #    _verbose        Verbosity
    # Available after _prepare
    #    _sourceF        FFT of source images
    #    _sourceR        real space source images
    #    _waveF          Current wavefunction (fourier space)
    #    _waveR          Current wavefunction (real space)
    #    _globalNorm     Ratio between series and factors
    def __init__(self, series, defocus=None):
        # Test param
        if len(series.shape) != 2:
            raise ValueError("Series must be 2D.")
        # Init Object
        self._series = series
        self._indexShape = series.indexShape
        self._indexSize = series.indexSize
        self._grid = series.grid
        self._initGrid()
        self._initDefocus(defocus)
        self._waveLength = calcWaveLength(series.attrs["voltage(kV)"])
        self._error = np.zeros(series.indexShape, dtype=float)
        self._shift = np.zeros(series.indexShape + (2,), dtype=float)
        self._tilt = np.zeros(series.indexShape + (2,), dtype=float)
        self._factor = np.ones(series.indexShape, dtype=complex)
        self._iteration = 0  # Current iteratiuon
        self._convergence = []
        self._lock = None
        # Public attributes
        self.verbose = 0
        self.adjustDefocus = True
        self.adjustShift = True
        self.adjustTilt = True
        self.margin = 1.0 / 6.0

    def _initGrid(self):
        qy, qx = self._grid.getRcprGrid()
        self._qy = qy
        self._qx = qx
        self._qx2 = qx ** 2
        self._qy2 = qy ** 2
        ry, rx = self._grid.getRealGrid()
        self._ry = ry
        self._rx = rx

    def _initDefocus(self, defocus):
        if defocus is None:
            defocus = np.zeros(self._indexShape, dtype=float)
        elif defocus.shape != self._indexShape:
            raise ValueError("defocus array must have same size as series.")
        self._defocus = defocus.copy()

    def _prepare(self):
        self._waveF = np.zeros(self._grid.shape, dtype=self._grid.complexType)
        self._sourceR = np.empty(self._indexShape + self._grid.shape, dtype=self._grid.complexType)
        self._grid.prepareFFT()
        space = self._series.attrs.get("space", +1)
        for n in range(self._indexSize):
            index = np.unravel_index(n, self._indexShape)
            # get FFT
            if space >= 0:
                self._sourceR[index, ...] = self._series[index].array
            else:
                self._sourceR[index, ...] = self._grid.backwardFFT(self._series[index].array)
            self._factor[index] = np.mean(self._sourceR[index])
        self._globalNorm = np.mean(abs(self._factor))
        self._factor /= self._globalNorm
        # print "Globalnorm:", self._globalNorm

    def _averageSingle(self, index):
        tmp = self._propagateSingle(self._sourceR[index], self._shift[index], self._defocus[index], self._tilt[index], self._factor[index])
        self._waveF += tmp

    def _evaluateSingleError(self, index):
        tmp = self._backPropagateSingle(self._waveF, self._shift[index], self._defocus[index], self._tilt[index], self._factor[index])
        self._error[index] = self._calcError(tmp, self._sourceR[index])

    def _dumpCurrentError(self, errorSum):
        if self.verbose > 1:
            T = self._grid.realSampling
            invT = np.linalg.inv(T)
            print("Optimizing after iteration %d" % self._iteration)
            print("\t[NN] sx[px] sy[px] tx[1/px] ty[1/px] def[nm] Ampl   Phase  Error")
            for n in range(self._indexSize):
                index = np.unravel_index(n, self._indexShape)
                shift_px = np.dot(invT, self._shift[index + (Ellipsis,)])
                tilt_px = np.dot(T, self._tilt[index + (Ellipsis,)])
                a = self._factor[index].real ** 2 + self._factor[index].imag ** 2
                p = atan2(self._factor[index].imag, self._factor[index].real)
                print("\t[%02d] %6.3f %6.3f %8.5f %8.5f %7.3f %6.4f %+6.3f %8e" % (n, shift_px[1], shift_px[0], tilt_px[1], tilt_px[0], self._defocus[index], a, p, self._error[index]))
        if self.verbose > 0:
            print("Iteration %3d: total error=%e" % (self._iteration, errorSum))

    def _propagateSingleFFT(self, sourceF, shift, defocus, factor):
        tmp = aligned_copy(sourceF)
        preD = pi * self._waveLength * defocus
        preSX = -2.0 * pi * shift[1]
        preSY = -2.0 * pi * shift[0]
        propX = np.exp(1.0j * (self._qx * preSX + self._qx2 * preD))
        propY = np.exp(1.0j * (self._qy * preSY + self._qy2 * preD))
        tmp *= propX * propY / factor
        return tmp

    def _backPropagateSingleFFT(self, sourceF, shift, defocus, factor):
        return self._propagateSingleFFT(sourceF, -shift, -defocus, 1.0 / factor)

    def _propagateSingle(self, sourceR, shift, defocus, tilt, factor):
        tmp = aligned_copy(sourceR)
        preTX = 2.0 * pi * tilt[1]
        preTY = 2.0 * pi * tilt[0]
        facX = np.exp(-1.0j * (self._rx * preTX))
        facY = np.exp(-1.0j * (self._ry * preTY))
        tmp *= facX * facY
        tmp = self._grid.forwardFFT(tmp)
        return self._propagateSingleFFT(tmp, shift, defocus, factor)

    def _backPropagateSingle(self, sourceF, shift, defocus, tilt, factor):
        tmp = self._backPropagateSingleFFT(sourceF, shift, defocus, factor)
        tmp = self._grid.backwardFFT(tmp)
        preTX = 2.0 * pi * tilt[1]
        preTY = 2.0 * pi * tilt[0]
        facX = np.exp(+1.0j * (self._rx * preTX))
        facY = np.exp(+1.0j * (self._ry * preTY))
        tmp *= facX * facY
        return tmp

    def _calcError(self, model, data, leastsq_compatible=False):
        rx = int(model.shape[1] * self.margin)
        ry = int(model.shape[0] * self.margin)
        sub = (slice(ry, -ry), slice(rx, -rx))
        res = model[sub] - data[sub]
        if leastsq_compatible:
            return abs(res).ravel()
        else:
            return np.sum(res.real ** 2 + res.imag ** 2)

    def _optimizeFunc(self, p, index, optimizeShift, optimizeTilt, optimizeDefocus):
        factor = complex(p[0], p[1])
        if optimizeShift:
            shift = p[2:4]
            pIter = 4
        else:
            shift = self._shift[index]
            pIter = 2
        if optimizeTilt:
            tilt = p[pIter:pIter + 2]
            pIter += 2
        else:
            tilt = self._tilt[index]
        if optimizeDefocus:
            defocus = p[pIter]
        else:
            defocus = self._defocus[index]
        if np.allclose(tilt, 0):
            model = self._backPropagateSingleFFT(self._waveF, shift, defocus, factor)
            model = self._grid.backwardFFT(model)
        else:
            model = self._backPropagateSingle(self._waveF, shift, defocus, tilt, factor)
        return self._calcError(model, self._sourceR[index], leastsq_compatible=True)

    def _optimizeSingle(self, index, optimizeShift, optimizeTilt, optimizeDefocus):
        p0 = (self._factor[index].real, self._factor[index].imag)
        if optimizeShift:
            p0 = p0 + tuple(self._shift[index])
        if optimizeTilt:
            p0 = p0 + tuple(self._tilt[index])
        if optimizeDefocus:
            p0 = p0 + (self._defocus[index],)
        pX, ier = leastsq(self._optimizeFunc, p0, args=(index, optimizeShift, optimizeTilt, optimizeDefocus), ftol=1e-4, epsfcn=1e-4)
        if ier not in [1, 2, 3, 4]:
            raise RuntimeError("Fit did not converge: index=%d" % index)
        self._factor[index] = complex(pX[0], pX[1])
        if optimizeShift:
            self._shift[index, ...] = pX[2:4]
            pIter = 4
        else:
            pIter = 2
        if optimizeTilt:
            self._tilt[index, ...] = pX[pIter:pIter + 2]
            pIter += 2
        if optimizeDefocus:
            self._defocus[index] = pX[pIter]

    def _averageFirst(self):
        self._waveF[...] = 0.0
        for n in range(self._indexSize):
            index = np.unravel_index(n, self._indexShape)
            if n > 0:
                self._optimizeSingle(index, self.adjustShift, False, False)
            self._factor[index] /= abs(self._factor[index])
            self._waveF *= n
            self._averageSingle(index)
            self._waveF /= n + 1
        self._waveR = self._grid.backwardFFT(self._waveF.copy())
        for n in range(self._indexSize):
            index = np.unravel_index(n, self._indexShape)
            self._evaluateSingleError(index)
        errorSum = np.sum(self._error)
        self._dumpCurrentError(errorSum)
        return errorSum

    def _averageCurrent(self):
        self._factor /= np.mean(abs(self._factor))  # Keep average factor at amplitude=1
        self._waveF[...] = 0.0
        for n in range(self._indexSize):
            index = np.unravel_index(n, self._indexShape)
            self._averageSingle(index)
        self._waveF /= self._indexSize
        self._waveR = self._grid.backwardFFT(self._waveF.copy())
        for n in range(self._indexSize):
            index = np.unravel_index(n, self._indexShape)
            self._evaluateSingleError(index)
        errorSum = np.sum(self._error)
        self._dumpCurrentError(errorSum)
        return errorSum

    def _showCurrent(self):
        import matplotlib.pyplot as plt
        tmp = self._waveR
        plt.subplot(121)
        plt.imshow(abs(tmp))
        plt.subplot(122)
        plt.imshow(np.arctan2(tmp.imag, tmp.real))
        plt.show()

    def average(self, iterations):
        """Average series for *iteration* rounds"""
        # Prepare
        if (self._iteration == 0):
            self._prepare()
            errorSum = self._averageFirst()
            self._convergence.append(errorSum)
        while iterations > 0:
            if self.verbose > 2:
                self._showCurrent()
            # Optimize ?
            for n in range(self._indexSize):
                index = np.unravel_index(n, self._indexShape)
                self._optimizeSingle(index, self.adjustShift, self.adjustTilt, (self._iteration >= 4) and self.adjustDefocus)
            # Book-keeping
            self._iteration += 1
            iterations -= 1
            errorSum = self._averageCurrent()
            self._convergence.append(errorSum)
        # Done

    def getCurrentAsDataSet(self):
        data = DataSet(self._waveR.shape, self._waveR.dtype)
        data.array[...] = self._waveR
        data.attrs.update(self._series.attrs)
        T = ScaleMatrix(2)
        T.set(self._grid.realSampling[::-1, ::-1])
        data.attrs["dim_offset"] = [0.0, 0.0]
        data.attrs["dim_scale"] = T.getCompact(minRank=1)
        data.attrs["dim_unit"] = ["nm", "nm"]
        data.attrs["space"] = +1
        data.attrs["shift(nm)"] = self._shift
        data.attrs["tilt(1/nm)"] = self._tilt
        data.attrs["defocus(nm)"] = self._defocus
        data.attrs["factor"] = self._factor
        data.attrs["error"] = self._error
        data.attrs["convergence"] = self._convergence
        return data

    def getCurrent(self):
        """The reconstructed (real space) wave."""
        return self._waveR.copy()

    def getVariance(self):
        var = np.zeros(self._grid.shape, self._grid.floatType)
        for n in range(self._indexSize):
            index = np.unravel_index(n, self._indexShape)
            delta = self._waveR - self._grid.backwardFFT(
                self._propagateSingle(self._sourceR[index], self._shift[index], self._defocus[index], self._tilt[index], self._factor[index]))
            var += delta.real ** 2 + delta.imag ** 2
        return var / (self._indexSize - 1)

    def getVarianceAsDataSet(self):
        data = DataSet(self._grid.shape, self._grid.floatType)
        data.array[...] = self.getVariance()
        data.attrs.update(self._series.attrs)
        T = ScaleMatrix(2)
        T.set(self._grid.realSampling[::-1, ::-1])
        data.attrs["dim_offset"] = [0.0, 0.0]
        data.attrs["dim_scale"] = T.getCompact(minRank=1)
        data.attrs["dim_unit"] = ["nm", "nm"]
        data.attrs["space"] = +1
        data.attrs["shift(nm)"] = self._shift
        data.attrs["tilt(1/nm)"] = self._tilt
        data.attrs["defocus(nm)"] = self._defocus
        data.attrs["factor"] = self._factor
        data.attrs["error"] = self._error
        data.attrs["convergence"] = self._convergence
        return data


def holoAverage(series, defocus=None, iterations=7, adjustTilt=True, adjustDefocus=True, adjustShift=True, margin=None, verbose=2, variance=False):
    averager = HoloAveraging(series, defocus=defocus)
    averager.verbose = verbose
    averager.adjustDefocus = adjustDefocus
    averager.adjustTilt = adjustTilt
    averager.adjustShift = adjustShift
    if margin is not None:
        averager.margin = margin
    averager.average(iterations)
    if variance:
        return averager.getCurrentAsDataSet(), averager.getVarianceAsDataSet()
    else:
        return averager.getCurrentAsDataSet()

