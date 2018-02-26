#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from setuptools import setup


long_description = """
Holoaverage is a Python script for the reconstruction and averaging of series of off-axis electron holograms, 
typically recorded in transmission electron holograms. The averaging is performed iteratively, such that instabilities 
of the microscope, like specimen and biprism drifts, can be tracked and corrected between consecutive exposures.

The averaging process is speed up, when also the `pyFFTW <http://hgomersall.github.com/pyFFTW/>`_ package is installed. 
However, is not a requirement for holoaverage and thus not automatically installed by ``pip``.

The source for holoaverage can be found on `GitHub <https://github.com/niermann/holoaverage>`_. The documentation can
be found on `ReadTheDocs <https://holoaverage.readthedocs.io>`_.
"""

setup(
    name='holoaverage',
    version='1.0.0',
    license='GPLv3+',
    description='Reconstruction and averaging of off-axis electron holograms as obtained by transmission electron microscopes.',
    long_description=long_description,
    author='Tore Niermann',
    author_email='tore.niermann@tu-berlin.de',
    url='https://github.com/niermann/holoaverage',
    packages=["holoaverage"],
    provides=["holoaverage"],
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    install_requires=['numpy', 'scipy', 'h5py'],
    extras_require={
        'rst': ['docutils>=0.11'],
        'fftw': ['pyfftw'],
    },
    entry_points={
        'console_scripts': [
            'holoaverage = holoaverage.main:main',
        ]
    },
    project_urls={
        "Documentation": "https://holoaverage.readthedocs.io",
        "Source Code": "https://github.com/niermann/holoaverage",
    }
)
