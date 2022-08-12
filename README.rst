About
=====

Holoaverage is a Python program for the reconstruction and averaging of series of off-axis electron holograms
recorded in a transmission electron microscope. The averaging is performed iteratively, such that instabilities of
the microscope, like specimen and biprism drifts, can be tracked and corrected between consecutive exposures.

The program is written and maintained by Tore Niermann (email: tore.niermann@tu-berlin.de).

Details on the usage can be found in the documentation. The documentation can be
found at

    https://holoaverage.readthedocs.io

Citation
--------

When you use the program in your research work, please cite the paper describing the details of the averaging method.
The details can be found in

        | T. Niermann and M. Lehmann
        | Averaging scheme for atomic resolution off-axis electron holograms
        | Micron 63 (2014) 28-34
        | doi: `10.1016/j.micron.2014.01.008 <http://dx.doi.org/10.1016/j.micron.2014.01.008>`_

The BibTeX entry for the paper is:

.. code-block:: bib

    @article{Niermann2014,
        title = "Averaging scheme for atomic resolution off-axis electron holograms",
        journal = "Micron",
        volume = "63",
        pages = "28 - 34",
        year = "2014",
        doi = "https://doi.org/10.1016/j.micron.2014.01.008",
        author = "T. Niermann and M. Lehmann",
        keywords = "Off-axis electron holography, High-resolution transmission electron microscopy, Iterative reconstruction"
    }


.. _sec-installation:

Installation
------------

The program is a Python program. Thus, a working Python 3.X distribution is required for running it (support for Python
2.X has been dropped since version 1.1.8). Beside the Python interpreter it requires the following
packages:

    * numpy (see `<www.numpy.org>`_)
    * scipy (see `<www.scipy.org>`_)
    * h5py (see `<www.hypy.org>`_)
    * pyFFTW is optional; speeds up the averaging (see `<https://pypi.python.org/pypi/pyFFTW>`_)

The package is tested with following Python versions:
    * Python 3.6, numpy 1.11, scipy 0.19, h5py 2.7, PyFFTW 0.10
    * Python 3.8, numpy 1.23, scipy 1.9, h5py 3.7, PyFFTW 0.13

The package can be most conveniently installed using the ``pip`` package manager. Make sure you have Python installed
and the Python interpreter is in your path. Go to the command line and execute:

.. code-block:: none

    python3 -m pip install --upgrade holoaverage

Holoaverage leverages the pyFFTW package for speed. If pyfftw can not be installed you can still use holoaverage
without problems. You can install pyFFTW by

.. code-block:: none

    python3 -m pip install --upgrade pyfftw

Up to date source versions can be found on the GitHub site: https://github.com/niermann/holoaverage

Bug reporting
-------------

When `reporting a bug <https://github.com/niermann/holoaverage/issues>`_ please include:

    * Your operating system name and version.
    * Any details about your local setup that might be helpful in troubleshooting.
    * Detailed steps to reproduce the bug.

License
-------

.. code-block:: none

    Holoaverage, program for reconstruction and averaging of electron holograms
    Copyright (C) 2018-2022 Tore Niermann

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
