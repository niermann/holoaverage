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
import math

from .grid import ScaleMatrix

__all__ = ('FilterFunction', )


class FilterFunction(object):
    """Define mask/filter.

    The mask can be isotropic or anisotropic. For anisotropic masks use max_q2 to specify the
    the normalized frequency.

    If max_q2 is used, g = sqrt(qT dot inv(max_q2) dot q) is the normalized frequency. Otherwise,
    if max_q is used, g = q / max_q is the normalized frequency. Either max_q or max_q2 must
    be used.

    The mask type parameter describes the kind of mask used. The parameter is either
    a string or a tuple with a string as first member. The string gives the type of
    the mask, further elements of the tuple give parameters for the mask. Supported values for
    type:

        "BUTTERWORTH": mask_type[1] gives order (n) of filter. Filter function is f(g)=1/(1 + g^2n)

        "GAUSSIAN": f(g)=exp(-g^2/2)

        "EDGE": f(g)=(g < 1)

        "none": f(q)=1
    """
    MASKTYPE_BUTTERWORTH = "BUTTERWORTH"
    MASKTYPE_GAUSSIAN = "GAUSSIAN"
    MASKTYPE_EDGE = "EDGE"
    MASKTYPE_NONE = "NONE"

    def __init__(self, max_q=None, max_q2=None, mask_type=None):
        self._max_q2 = ScaleMatrix(2, True)
        if max_q2 is not None:
            if max_q:
                raise ValueError("Either max_q or max_q2 must be set.")
            self.max_q2 = max_q2
        elif max_q is not None:
            self.max_q = max_q
        else:
            raise ValueError("Either max_q or max_q2 must be set.")
        if mask_type is None:
            mask_type = self.MASKTYPE_EDGE
        self.mask_type = mask_type

    @property
    def max_q(self):
        """Cutoff frequency for isotropic filters in 1/nm. For anisotropic filters this is the
        square root of the average eigen value of max_q2"""
        return math.sqrt(self._max_q2.getEucledianMean())

    @max_q.setter
    def max_q(self, value):
        self._max_q2.set(value ** 2)

    @property
    def max_q2(self):
        """Cutoff frequency for anisotropic filters in 1/nm^2."""
        return self._max_q2.getMatrix()

    @max_q2.setter
    def max_q2(self, value):
        self._max_q2.set(value)

    @property
    def mask_type(self):
        """The mask type parameter."""
        return self._mask_type

    @mask_type.setter
    def mask_type(self, value):
        if isinstance(value, str):
            value = (value,)
        value = (value[0].upper(),) + tuple(value[1:])
        if value[0] == self.MASKTYPE_BUTTERWORTH:
            if len(value) != 2 or int(value[1]) < 0:
                raise ValueError("Butterworth filter requires non-negative order parameter")
        elif value[0] not in [self.MASKTYPE_GAUSSIAN, self.MASKTYPE_EDGE, self.MASKTYPE_NONE]:
            raise ValueError("Unsupported mask type: %s" % value[0])
        self._mask_type = value

    def calculate(self, grid, roi=None, origin=None):
        """Calculates mask/filter for given grid.

        The roi parameter allows to only create the filter for a subarea of the shape given by
        grid. The roi parameter refers to the unshifted grid!

        The origin parameter allows to give the (x,y) center of the filter explicitly (in 1/nm)
        """
        if len(grid.shape) != 2:
            raise ValueError("Grid must be 2D.")
        if roi is None:
            roi = [0, 0, grid.shape[1], grid.shape[0]]
        # Calculate normalized frequency
        qy, qx = grid.getRcprGrid()
        if origin is not None:
            qy -= float(origin[1])
            qx -= float(origin[0])
        if self._max_q2.isDiagonal():
            qmax2 = self._max_q2.getMatrix()
            g2 = qy ** 2 / qmax2[1, 1] + qx ** 2 / qmax2[0, 0]
        else:
            iqmax2 = np.linalg.inv(self._max_q2.getMatrix())
            g2 = qy ** 2 * iqmax2[1, 1] + qx ** 2 * iqmax2[0, 0] + 2 * qx * qy * iqmax2[0, 1]
        # Calculate filter function
        out = np.empty((roi[3] - roi[1], roi[2] - roi[0]), dtype=grid.floatType)
        if self._mask_type[0] == self.MASKTYPE_BUTTERWORTH:
            # Not real butterworth, but gain of butterworth
            # See http://de.wikipedia.org/wiki/Butterworth-Filter
            order = int(self._mask_type[1])
            denom = 1.0 + np.power(g2[roi[1]:roi[3], roi[0]:roi[2]], order)
            out = np.divide(1.0, denom, out=out)
        elif self._mask_type[0] == self.MASKTYPE_GAUSSIAN:
            out = np.exp(-0.5 * g2[roi[1]:roi[3], roi[0]:roi[2]], out)
        elif self._mask_type[0] == self.MASKTYPE_EDGE:
            out[...] = g2[roi[1]:roi[3], roi[0]:roi[2]] < 1.0
        else:
            out[...] = 1.0
        return out

    def calculate_shifted(self, grid, roi=None, origin=None):
        """Calculates mask/filter for given grid.

        The roi parameter allows to only create the filter for a subarea of the shape given by
        grid. The roi parameter refers to the fft-shifted grid!

        The origin parameter allows to give the (x,y) center of the filter explicitly (in 1/nm)
        """
        if len(grid.shape) != 2:
            raise ValueError("Grid must be 2D.")
        if roi is None:
            roi = [0, 0, grid.shape[1], grid.shape[0]]
        # Calculate normalized frequency
        qy, qx = grid.getRcprGrid()
        if origin is not None:
            qy -= float(origin[1])
            qx -= float(origin[0])
        if self._max_q2.isDiagonal():
            qmax2 = self._max_q2.getMatrix()
            g2 = np.fft.fftshift(qy**2 / qmax2[1, 1] + qx**2 / qmax2[0, 0])
        else:
            iqmax2 = np.linalg.inv(self._max_q2.getMatrix())
            g2 = np.fft.fftshift(qy**2 * iqmax2[1, 1] + qx**2 * iqmax2[0, 0] + 2 * qx * qy * iqmax2[0, 1])
        # Calculate filter function
        out = np.empty((roi[3] - roi[1], roi[2] - roi[0]), dtype=grid.floatType)
        if self._mask_type[0] == self.MASKTYPE_BUTTERWORTH:
            # Not real butterworth, but gain of butterworth
            # See http://de.wikipedia.org/wiki/Butterworth-Filter
            order = int(self._mask_type[1])
            denom = 1.0 + np.power(g2[roi[1]:roi[3], roi[0]:roi[2]], order)
            out = np.divide(1.0, denom, out=out)
        elif self._mask_type[0] == self.MASKTYPE_GAUSSIAN:
            out = np.exp(-0.5 * g2[roi[1]:roi[3], roi[0]:roi[2]], out)
        elif self._mask_type[0] == self.MASKTYPE_EDGE:
            out[...] = g2[roi[1]:roi[3], roi[0]:roi[2]] < 1.0
        else:
            out[...] = 1.0
        return out
