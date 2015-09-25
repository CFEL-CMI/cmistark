#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# Copyright (C) 2008,2009,2012,2013,2014,2015 Jochen Küpper <jochen.kuepper@cfel.de>


import numpy,os
from setuptools import setup
from Cython.Build import cythonize

extra_compile_args = []
library_dirs = []

long_description = """CMI Python extensions for Stark effect calculations

Python extensions for Stark effect calculations of molecules.
Developed by Jochen Küpper and the Controlled Molecule Imaging Group at the Center for Free-Electron Laser Science,
Hamburg, Germany.

Original author:    Jochen Küpper <jochen.kuepper@cfel.de>
Current maintainer: Jochen Küpper <jochen.kuepper@cfel.de>
See the distribution files AUTHORS and THANKS for further contributions.
"""


setup(name="cmistark",
      author              = "Jochen Küpper, CFEL-CMI group, et al (see AUTHORS)",
      author_email        = "jochen.kuepper@cfel.de",
      maintainer          = "Jochen Küpper and the CFEL-CMI group",
      maintainer_email    = "jochen.kuepper@cfel.de",
      url                 = "http://desy.cfel.de/cid/cmi/cmistark",
      description         = "CMI Python Stark effect extensions",
      version             = "1.1.dev0",
      long_description    = long_description,
      license             = "GPL",
      ext_modules         = cythonize('cmistark/*.pyx'), include_dirs=[numpy.get_include()],
      scripts             = ['scripts/cmistark_brute-force-orientation',
                             'scripts/cmistark_calculate_energy',
                             'scripts/cmistark_plot_energy',
                             'scripts/cmistark_print_energy'],
      test_suite          = 'tests',
      )
