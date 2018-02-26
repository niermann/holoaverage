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

import numpy as np

__all__ = ('pyfftw_present', 'zeros_aligned', 'empty_aligned', 'ones_aligned', 'aligned_copy',
           'fromHermite', 'create_transform_wrapper')


_pyfftw_present = None
zeros_aligned = None
empty_aligned = None
ones_aligned = None


def disable_pyfftw():
    # Enable fallbacks
    global _pyfftw_present, zeros_aligned, empty_aligned, ones_aligned

    _pyfftw_present = False
    zeros_aligned = np.zeros
    empty_aligned = np.empty
    ones_aligned = np.ones


def enable_pyfftw():
    # Use aligned allocators
    global _pyfftw_present, zeros_aligned, empty_aligned, ones_aligned

    import pyfftw
    _pyfftw_present = True
    zeros_aligned = pyfftw.zeros_aligned
    empty_aligned = pyfftw.empty_aligned
    ones_aligned  = pyfftw.ones_aligned


try:
    enable_pyfftw()
except ImportError:
    disable_pyfftw()


def pyfftw_present():
    """Return whether pyfftw is present"""
    global _pyfftw_present
    return _pyfftw_present


def aligned_copy(a):
    """
    Create an aligned copy of the input array.

    :param a: Input array.
    :returns: Aligned array
    """
    result = empty_aligned(a.shape, dtype=a.dtype)
    result[...] = a
    return result


def fromHermite(arr, shape, axis=-1):
    '''
    Returns full hermitian array from lower half.

    The ``rfft``-like FFTs will return only one half of the transformed
    array, namely the positive frequencies (or frequencies up to Nyquist), since
    the other half can be calculated from the hermitian nature of the real data FT.

    This function will expand such an hermitian half array into the full array. The shape
    of the full array is required, since the shape of the half-array is ambiguous.

    :param arr: Array to expand
    :param shape: Shape of expanded array
    :param axis: On which axis the array is to expanded (defaults to last)

    :returns: Expanded array
    '''
    arr = np.atleast_1d(arr)
    if axis < 0:
        axis = arr.ndim + axis
        if axis < 0:
            raise ValueError('Invalid axis.')
    shape = np.atleast_1d(shape)

    # Test shapes
    if arr.ndim != len(shape):
        raise ValueError('Number of dimensions between input array and output shape do not match.')
    for n in range(arr.ndim):
        if n != axis:
            if shape[n] != arr.shape[n]:
                raise ValueError("Shape of input and output arrays do not match.")
        else:
            if (shape[axis] // 2) + 1 != arr.shape[axis]:
                raise ValueError("Array 'arr' cannot be expanded to demanded shape.")

                # Copy "left" half
    output = np.empty(shape, dtype=arr.dtype)
    index = [slice(None)] * arr.ndim
    index[axis] = slice(0, shape[axis] // 2 + 1)
    output[index] = arr

    # "right" half is conjugated
    oIndex = [slice(None)] * arr.ndim
    iIndex = oIndex[:]
    oIndex[axis] = slice(shape[axis] // 2 + 1, None)
    iIndex[axis] = slice((shape[axis] - 1) // 2, 0, -1)
    output[oIndex] = np.conj(arr[iIndex])

    # Swap halfs on ride side on other directions
    for n in range(arr.ndim):
        if n == axis:
            continue
        index = [0] + list(range(shape[n] - 1, 0, -1))
        output[oIndex] = np.take(output[oIndex], index, n)
    return output


def create_transform_wrapper(transform):
    """Return wrapper for Fourier transformations with output argument"""
    def inner(x, out):
        tmp = transform(x)
        np.copyto(out, tmp)
        return tmp
    return inner