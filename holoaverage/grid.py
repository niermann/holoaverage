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

from .fft import empty_aligned, pyfftw_present, create_transform_wrapper

__all__ = ("ScaleMatrix", "Grid")


EPSILON = 1e-6


class ScaleMatrix(object):
    """
    Matrix that can be initialized by a scalar (results in diagonal with scalar value),
    by a vector (the diagonal), or by a matrix
    """
    __slots__ = ['_ndim', '_symmetric', '_value', '_matrix']

    def __init__(self, ndim, symmetric=False):
        self._ndim = ndim
        self._symmetric = symmetric
        self.set(None)

    @property
    def ndim(self):
        """Number of dimensions."""
        return self._ndim

    def set(self, value):
        """Set the value. value can be either a scalar, a vector or a matrix. Returns self."""
        if value is None:
            self._value = None
            self._matrix = None
            return self
        if isinstance(value, ScaleMatrix):
            value = value.getMatrix()
        else:
            value = np.array(value, dtype=float, copy=True)
        if value.size == 1:
            self._value = np.asscalar(value)
            self._matrix = np.eye(self._ndim) * self._value
        elif value.ndim == 1:
            if (value.size != self._ndim):
                raise ValueError("Invalid shape, expected (%d) or (%d,%d)." % ((self._ndim,) * 3))
            self._value = value
            self._matrix = np.diag(value)
        else:
            if (value.shape != (self._ndim,) * 2):
                raise ValueError("Invalid shape, expected (%d) or (%d,%d)." % ((self._ndim,) * 3))
            if self._symmetric and not np.allclose(value.T, value):
                raise ValueError("Matrix must be symmetric.")
            self._value = value
            self._matrix = value
        return self

    def setReverse(self, value):
        """Set the value. value can be either a scalar, a vector or a matrix,
        the order of indices is reversed. Returns self."""
        if value is None:
            self._value = None
            self._matrix = None
            return self
        if isinstance(value, ScaleMatrix):
            value = value.getMatrix()
        else:
            value = np.array(value, dtype=float, copy=True)
        if value.size == 1:
            self._value = np.asscalar(value)
            self._matrix = np.eye(self._ndim) * self._value
        elif value.ndim == 1:
            if (value.size != self._ndim):
                raise ValueError("Invalid shape, expected (%d) or (%d,%d)." % ((self._ndim,) * 3))
            self._value = value
            self._matrix = np.diag(value[::-1])
        else:
            if (value.shape != (self._ndim,) * 2):
                raise ValueError("Invalid shape, expected (%d) or (%d,%d)." % ((self._ndim,) * 3))
            if self._symmetric and not np.allclose(value.T, value):
                raise ValueError("Matrix must be symmetric.")
            self._value = value
            self._matrix = value[::-1, ::-1]
        return self

    def get(self):
        """Returns scaling."""
        return self._value

    def getReverse(self, obj):
        """Returns scaling with reversed index order."""
        if not self._value:
            return None
        if np.isscalar(self._value):
            return self._value
        elif self._value.ndim == 1:
            return self._value[::-1]
        else:
            return self._value[::-1, ::-1]

    def getCompact(self, minRank=0):
        """Returns scaling in most compact (scalar, vector, matrix) form.
        minRank allows to adjust the most compact form returned (0=scalar, 1=vector, 2=matrix)"""
        A = self._matrix
        if A is None:
            return None
        D = np.diag(A)
        if minRank >= 2 or not np.allclose(np.diag(D), A):
            return A
        if minRank >= 1 or not np.allclose(D, D[0]):
            return D
        else:
            return float(D[0])

    def getMatrix(self):
        """Returns scaling in matrix form."""
        A = self._matrix
        if A is None:
            return None
        else:
            return A.copy()

    def isSymmetric(self):
        A = self._matrix
        if A is None:
            return False
        return np.allclose(A.T, A)

    def isDiagonal(self):
        A = self._matrix
        if A is None:
            return False
        return np.allclose(np.diag(np.diag(A)), A)

    def isUniform(self):
        A = self._matrix
        if A is None:
            return False
        D = np.diag(A)
        if not np.allclose(np.diag(D), A):
            return False
        return np.allclose(D, D[0])

    def getGeometricMean(self):
        A = self._matrix
        if A is None:
            return 0.0
        return np.power(np.linalg.det(A), 1.0 / float(self._ndim))

    def getEucledianMean(self):
        A = self._matrix
        if A is None:
            return 0.0
        if not self._symmetric:
            raise ValueError("What are you doing? You might want to use the geometric mean")
        return np.mean(np.linalg.eigvalsh(A))


class Grid(object):
    """Helper class for grided calculations."""

    # Private members
    #   _realSampling    np.ndarray[ndim,ndim]: Real space sampling in for 2D: ((YY,YX),(XY,XX)) in nm
    #   _rcprSampling    np.ndarray[ndim,ndim]: Reciprocal space sampling in 2D: ((YY,YX),(XY,XX)) in 1/nm
    #   _shape           tuple: Size of simulation grid, len=ndim
    #   _complexType     np.dtype: Complex type
    #   _floatType       np.dtype: Float type
    #   _planFFT         FFT plan (destroys input)
    #   _planIFFT        IFFT plan (destroys input)
    #   _fft_nthreads    Number of threads to use in FFT plans
    def __init__(self, shape, realSampling=None, rcprSampling=None, dtype=np.complex64):
        """
        Initializes grid. Either realSampling or rcprSampling must be given.

        Arguments:
            shape
                Tuple of int: Shape of grid
            realSampling
                Float: Isotropic sampling in nm
                1D-array: Diagonal sampling in nm
                2D-array: Anisotropic sampling in nm
            rcprSampling
                Float: Isotropic sampling in 1/nm
                1D-array: Diagonal sampling in 1/nm
                2D-array: Anisotropic sampling 1/nm
            dtype
                np.dtype: Data type to use
        """
        self._shape = tuple(shape)
        dtype = np.dtype(dtype)
        if dtype.kind == 'c':
            self._complexType = np.dtype('c%d' % dtype.itemsize)
            self._floatType = np.dtype('f%d' % (dtype.itemsize / 2))
        elif dtype.kind == 'f':
            self._complexType = np.dtype('c%d' % (dtype.itemsize * 2))
            self._floatType = np.dtype('f%d' % dtype.itemsize)
        else:
            raise TypeError("Unsupported data type.")
        self._realSampling = ScaleMatrix(2)
        self._rcprSampling = ScaleMatrix(2)
        if realSampling is not None:
            self._setRealSampling(realSampling)
        elif rcprSampling is not None:
            self._setRcprSampling(rcprSampling)
        else:
            raise ValueError("Either 'realSampling' or 'rcprSampling' must be given.")
        self._planFFT = None
        self._planIFFT = None
        self._planRFFT = None
        self._planIRFFT = None

    def __getstate__(self):
        return {'shape': self.shape,
                'complexType': self._complexType,
                'floatType': self._floatType,
                'realSampling': self._realSampling,
                'rcprSampling': self._rcprSampling}

    def __setstate__(self, state):
        self._shape = state["shape"]
        self._complexType = state["complexType"]
        self._floatType = state["floatType"]
        self._realSampling = state["realSampling"]
        self._rcprSampling = state["rcprSampling"]
        self._planFFT = None
        self._planIFFT = None
        self._planRFFT = None
        self._planIRFFT = None

    @staticmethod
    def fromDataSet(dataset):
        """Creates grid object from DataSet."""
        shape = dataset.shape
        dtype = dataset.dtype
        if dtype.kind in ['i', 'u']:
            dtype = np.dtype('f%d' % max(dtype.itemsize, 4))
        elif dtype.kind not in ['f', 'c']:
            raise ValueError("Unsupported data type.")
        if len(dataset.shape) != 2:
            raise ValueError("Expected 2D dataset.")

        if 'dim_unit' in dataset.attrs:
            dim_unit = dataset.attrs['dim_unit'][::-1]
            if dim_unit[0] != dim_unit[1]:
                raise ValueError("Expected uniform units in both dimensions.")
        else:
            dim_unit = None
        dim_scale = ScaleMatrix(len(dataset.shape))
        if 'dim_scale' in dataset.attrs:
            dim_scale.set(dataset.attrs['dim_scale']).getMatrix()
        else:
            dim_scale.set(np.ones(len(dataset.shape)))

        space = dataset.attrs.get('space', 0)
        if dim_unit is not None:
            if space == 0:
                if dim_unit[0] == "1/nm":
                    space = -1
                else:
                    space = +1
            elif space > 0:
                if dim_unit[0] != 'nm':
                    warnings.warn("Expected unit to be 'nm' not '%s'" % dim_unit[0])
            else:
                if dim_unit[0] != '1/nm':
                    warnings.warn("Expected unit to be '1/nm' not '%s'" % dim_unit[0])

        if space >= 0:
            return Grid(shape, realSampling=dim_scale, dtype=dtype)
        else:
            dim_scale.set(dim_scale.getMatrix() * np.atleast_1d(shape))
            return Grid(shape, rcprSampling=dim_scale, dtype=dtype)

    @property
    def complexType(self):
        """Complex type of grid."""
        return self._complexType

    @property
    def floatType(self):
        """Float type of grid."""
        return self._floatType

    @property
    def realSampling(self):
        """Real space sampling, as matrix: index [i,j] corresponds to position dot(realSampling, [i,j])."""
        return self._realSampling.getMatrix()

    def _setRealSampling(self, value):
        self._realSampling.set(value)
        self._rcprSampling.set(np.linalg.inv(self._realSampling.getMatrix()))

    @property
    def rcprSampling(self):
        """Reciprocal space sampling, as matrix: index [i,j] corresponds to position dot(rcprSampling, [i,j]/shape)."""
        return self._rcprSampling.getMatrix()

    def _setRcprSampling(self, value):
        self._rcprSampling.set(value)
        self._realSampling.set(np.linalg.inv(self._rcprSampling.getMatrix()))

    def hasDiagonalSampling(self):
        """Returns whether the sampling matrix is diagonal (sampling along grid axes)."""
        return self._realSampling.isDiagonal()

    def hasUniformSampling(self):
        """Returns whether the sampling is isotropic (same along grid axes)."""
        return self._realSampling.isUniform()

    @property
    def shape(self):
        """Shape of grid."""
        return self._shape

    @property
    def hermiteShape(self):
        """Shape of hermite grid."""
        return self._shape[:-1] + (self._shape[-1] // 2 + 1,)

    @property
    def ndim(self):
        """Rank of grid (number of dimensions)."""
        return len(self._shape)

    @property
    def size(self):
        """Total number of elements in grid."""
        return np.prod(self._shape)

    def prepareFFT(self, forwardFFT=True, backwardFFT=True, forwardRFFT=False, backwardRFFT=False, verbose=0, pyfftw_planner="FFTW_MEASURE", pyfftw_threads=1):
        """Prepare FFTs by doing any precalculation needed for them.

        Depending on the PyFFT backend, this step might need a few seconds. If the FFTs are not prepared by a call
        to prepareFFT, they are lazily prepared when needed. You can fine-tune, which FFTs are actually needed using the
        boolean flags.

        :param forwardFFT: Prepare forward FFT
        :param backwardFFT: Prepare backward FFT
        :param forwardRFFT: Prepare (forward) real to hermite FFT
        :param backwardFFT: Prepare (backward) hermite to real FFT
        :param verbose: Set to >0, if preparation time should be displayed.
        :param pyfftw_planner: Planner flags to use for pyfftw
        :param pyfftw_threads: Number of threads to use for pyfftw
        """
        if pyfftw_present():
            import pyfftw
            if verbose > 0:
                print("Planning FFT (FFT=%s, IFFT=%s, RFFT=%s, IRFFT=%s) ..." % (forwardFFT, backwardFFT, forwardRFFT, backwardRFFT))

            from time import clock
            timer = clock()

            if forwardFFT and self._planFFT is None:
                tmp = empty_aligned(self.shape, self.complexType)
                self._planFFT = pyfftw.builders.fft2(tmp, overwrite_input=True, planner_effort=pyfftw_planner, threads=pyfftw_threads, avoid_copy=True)
            if backwardFFT and self._planIFFT is None:
                tmp = empty_aligned(self._shape, dtype=self.complexType)
                self._planIFFT = pyfftw.builders.ifft2(tmp, overwrite_input=True, planner_effort=pyfftw_planner, threads=pyfftw_threads, avoid_copy=True)
            if forwardRFFT and self._planRFFT is None:
                tmp = empty_aligned(self._shape, dtype=self.floatType)
                self._planRFFT = pyfftw.builders.rfft2(tmp, overwrite_input=True, planner_effort=pyfftw_planner, threads=pyfftw_threads, avoid_copy=True)
            if backwardRFFT and self._planIRFFT is None:
                tmp = empty_aligned(self.hermiteShape, dtype=self.complexType)
                self._planIRFFT = pyfftw.builders.irfft2(tmp, self.shape, overwrite_input=True, planner_effort=pyfftw_planner, threads=pyfftw_threads, avoid_copy=True)

            timer = clock() - timer
            if verbose > 0:
                print("    Took %f seconds." % timer)
        else:
            self._planFFT   = create_transform_wrapper(np.fft.fft2)
            self._planIFFT  = create_transform_wrapper(np.fft.ifft2)
            self._planRFFT  = create_transform_wrapper(np.fft.rfft2)
            self._planIRFFT = create_transform_wrapper(np.fft.irfft2)

    def forwardFFT(self, inputArray, out=None):
        """Performs forward FFT (might destroy inputArray)"""
        if self._planFFT is None:
            self.prepareFFT(forwardFFT=True, backwardFFT=False, forwardRFFT=False, backwardRFFT=False)
        if out is None:
            out = empty_aligned(self.shape, self.complexType)
        self._planFFT(inputArray, out)
        return out

    def backwardFFT(self, inputArray, out=None):
        """Performs backward FFT (might destroy inputArray)"""
        if self._planIFFT is None:
            self.prepareFFT(forwardFFT=False, backwardFFT=True, forwardRFFT=False, backwardRFFT=False)
        if out is None:
            out = empty_aligned(self.shape, self.complexType)
        self._planIFFT(inputArray, out)
        return out

    def forwardRFFT(self, inputArray, out=None):
        """Performs forward FFFT (might destroy inputArray)"""
        if self._planRFFT is None:
            self.prepareFFT(forwardFFT=False, backwardFFT=False, forwardRFFT=True, backwardRFFT=False)
        if out is None:
            out_shape = self.shape[:-1] + (self.shape[-1] // 2 + 1,)
            out = empty_aligned(out_shape, self.complexType)
        self._planRFFT(inputArray, out)
        return out

    def backwardRFFT(self, inputArray, out=None):
        """Performs backward FFFT (might destroy inputArray)"""
        if self._planIRFFT is None:
            self.prepareFFT(forwardFFT=False, backwardFFT=False, forwardRFFT=False, backwardRFFT=True)
        if out is None:
            out = empty_aligned(self.shape, self.floatType)
        self._planIRFFT(inputArray, out)
        return out

    def getHermiteFreq(self):
        """Returns hermite grid of pixel frequencies, for current shape"""
        result = []
        ndim = self.ndim
        for n, dim in enumerate(self._shape[:-1]):
            idx = np.arange(dim, dtype=self.floatType)
            idx[0:(dim + 1) // 2] = np.arange((dim + 1) // 2, dtype=self.floatType)
            idx[-dim // 2:] = np.arange(-dim // 2, 0, dtype=self.floatType)
            idx /= float(dim)
            s = (1,) * n + (dim,) + (1,) * (ndim - n - 1)
            result.append(idx.reshape(s))
        idxN = np.arange(self._shape[-1] // 2 + 1, dtype=self.floatType) / float(self._shape[-1])
        result.append(idxN)
        return result

    def getFreq(self):
        """Returns full grid of pixel frequencies, for current shape"""
        result = []
        ndim = self.ndim
        for n, dim in enumerate(self._shape):
            idx = np.arange(dim, dtype=self.floatType)
            idx[0:(dim + 1) // 2] = np.arange((dim + 1) // 2, dtype=self.floatType)
            idx[-dim // 2:] = np.arange(-dim // 2, 0, dtype=self.floatType)
            idx /= float(dim)
            s = (1,) * n + (dim,) + (1,) * (ndim - n - 1)
            result.append(idx.reshape(s))
        return result

    def getRcprGrid(self):
        """Returns full grid of reciprocal frequencies, for current shape and sampling"""
        idx = self.getFreq()
        T = self.rcprSampling
        if self.hasDiagonalSampling():
            result = idx
            Td = np.diag(T)
            for i in range(self.ndim):
                result[i] *= Td[i]
        else:
            result = np.zeros((self.ndim,) + self._shape, dtype=self.floatType)
            for i in range(self.ndim):
                for j in range(self.ndim):
                    result[i] = result[i] + T[i, j] * idx[j]
        return result

    def getRealGrid(self):
        """Returns full grid of real coordinates, for current shape and sampling"""
        idx = []
        ndim = self.ndim
        for n, dim in enumerate(self._shape):
            tmp = np.arange(dim, dtype=self.floatType)
            s = (1,) * n + (dim,) + (1,) * (ndim - n - 1)
            idx.append(tmp.reshape(s))

        T = self.realSampling
        if self.hasDiagonalSampling():
            result = idx
            Td = np.diag(T)
            for i in range(ndim):
                result[i] *= Td[i]
        else:
            result = np.zeros((ndim,) + self._shape, dtype=self.floatType)
            for i in range(ndim):
                for j in range(ndim):
                    result[i] = result[i] + T[i, j] * idx[j]
        return result
