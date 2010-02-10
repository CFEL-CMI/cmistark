#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008,2009 Jochen Küpper <software@jochen-kuepper.de>

from numpy.distutils.core import setup, Extension
import sys
import os

extra_compile_args = []
library_dirs = []

long_description = """JK Python extensions for Stark effect calculations

Python extensions for Stark effect calculations of molecules.

Original authors:   Jochen Küpper <software@jochen-kuepper.de>
Current maintainer: Jochen Küpper <software@jochen-kuepper.de>
"""

version_major = sys.version_info[0]
version_minor = sys.version_info[1]

setup(name="jkstark",
      author              = "Jochen Küpper",
      author_email        = "software@jochen-kuepper.de",
      description         = "JK Python Stark effect extensions",
      license             = "GPL",
      url                 = "http://python.jochen-kuepper.de",
      version             = "0.0.1",
      long_description    = long_description,
      package_dir         = {'jkstark': 'lib'},
      packages            = ['jkstark'],
      ext_modules         = [Extension('jkstark._wigner_gsl',
                                       sources = ['src/wigner_gsl.c', 'src/wigner_gsl_module.c'],
                                       libraries = ['gsl', 'gslcblas']),
                             Extension('jkstark._wigner_avda',
                                       sources = ['src/wigner_avda.f',],
                                       extra_link_args = ['-bundle'],
                                       libraries = ['dl', 'python%d.%d' % (version_major,
                                                                           version_minor)]),
                             Extension('jkstark._wigner_fft',
                                       sources = ['src/wigner_fft.f'],
                                       include_dirs = ['/opt/local/include', '/usr/local/include', '/usr/include'],
                                       extra_link_args = ['-bundle'],
                                       libraries = ['fftw3f', 'fftw3',
                                                    'dl', 'python%d.%d' % (version_major,
                                                                           version_minor)]),
                             ],
      scripts             = ['scripts/jkstark_brute-force-orientation',
                             'scripts/jkstark_alignment_orientation',
                             'scripts/jkstark_calculate_energy',
                             'scripts/jkstark_plot_energy',
                             'scripts/jkstark_plot_energy_linear',
                             'scripts/jkstark_print_energy',
                             'scripts/plot_symtopwavefunction.py',
                             'scripts/plot_asymtopwavefunction.py']
      )

