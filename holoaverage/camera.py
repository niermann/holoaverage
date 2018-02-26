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

__ALL__ = ['ParameterizedMTF']


class ParameterizedMTF(object):
    """
    MTF parameterization.

    The parameterized MTF :math:`M(q)` consists of a list of terms :math:`f_n(q)`, which
    are summed up to result in the MTF:

    .. math::
        M(q) = \\sum_n f_n(q)

    The terms itself are given as tuples, where the first element of the tuple always
    identifies the type of the term. q is always given in 1/px, Nyquist is at 0.5 1/px

    Allowed values are:

        * (ParameterizedMTF.TERM_CONSTANT, A)

        .. math::
            f(q) = A

        * (ParameterizedMTF.TERM_GAUSSIAN, A, B):

        .. math::
            f(q) = A \\exp(-B q^2)

        * (ParameterizedMTF.TERM_LORENTZIAN, A, B):

        .. math::
            f(q) = A / (B + q^2)

    The :attr:`mtf_param` is a lists of these terms.
    """

    TERM_CONSTANT = "CONSTANT"
    TERM_GAUSSIAN = "GAUSSIAN"
    TERM_LORENTZIAN = "LORENTZIAN"

    def __init__(self, mtf_param):
        # Check MTF param
        self.mtf_param = []
        for term in mtf_param:
            term_type = term[0].upper()
            if term_type == ParameterizedMTF.TERM_CONSTANT:
                self.mtf_param.append((ParameterizedMTF.TERM_CONSTANT, float(term[1])))
            elif term_type == ParameterizedMTF.TERM_GAUSSIAN:
                self.mtf_param.append((ParameterizedMTF.TERM_GAUSSIAN, float(term[1]), float(term[2])))
            elif term_type == ParameterizedMTF.TERM_LORENTZIAN:
                self.mtf_param.append((ParameterizedMTF.TERM_LORENTZIAN, float(term[1]), float(term[2])))
            else:
                raise ValueError("Unexpected constant in MTF parameterization.")

    def mtf(self, shape, binning=1, scale=1.0, dtype=np.float32):
        """
        Returns the MTF.

        :param shape: Shape of the array to return (corresponds to size of ROI)
        :type shape: tuple of int
        :param binning: Binning (if tuple for X,Y direction)
        :type binning: int / tuple of int
        :param scale: Scaling (relative to binned pixels, i.e. 2 means two mtf-pixels are one binned pixel)
        :type scale: int / tuple of int
        :param dtype: Type of output array
        :type dtype: np.dtype
        """
        if np.isscalar(binning):
            binning = np.array((binning, binning), dtype=float)
        else:
            binning = np.array(binning, dtype=float)
        if np.isscalar(shape):
            shape = (shape, shape)
        if np.isscalar(scale):
            scale = np.array((scale, scale), dtype=float)
        else:
            scale = np.array(scale, dtype=float)
        qx = np.fft.fftfreq(shape[1], binning[0] / scale[0]).astype(dtype, copy=False).reshape(1, -1)
        qy = np.fft.fftfreq(shape[0], binning[1] / scale[1]).astype(dtype, copy=False).reshape(-1, 1)
        q2 = qx**2 + qy**2
        result = np.zeros(shape, dtype=dtype)
        for term in self.mtf_param:
            if term[0] == ParameterizedMTF.TERM_CONSTANT:
                result += term[1]
            elif term[0] == ParameterizedMTF.TERM_GAUSSIAN:
                result += term[1] * np.exp(-term[2] * q2)
            elif term[0] == ParameterizedMTF.TERM_LORENTZIAN:
                result += term[1] / (term[2] + q2)
            else:
                raise ValueError("Unknown parameter %s" % repr(term))
        result *= np.sinc(qx * binning[0]) * np.sinc(qy * binning[1])
        return result
