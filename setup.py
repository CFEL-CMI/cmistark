#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# Copyright (C) 2008,2009,2012,2013,2014 Jochen Küpper <jochen.kuepper@cfel.de>


import os
from setuptools import setup

extra_compile_args = []
library_dirs = []

long_description = """CMI Python extensions for Stark effect calculations

Python extensions for Stark effect calculations of molecules.
Developed by the Controlled Molecule Imaging Group at the Center for Free-Electron Laser Science, Hamburg, Germany.

Original author:    Jochen Küpper <jochen.kuepper@cfel.de>
Current maintainer: Yuan-Pin Chang <yuan.pin.chang@cfel.de>
See the distribution files AUTHORS and THANKS for further contributions.
"""


setup(name="cmistark",
      author              = "Jochen Küpper, CFEL-CMI group, et al (see AUTHORS)",
      author_email        = "jochen.kuepper@cfel.de",
      description         = "CMI Python Stark effect extensions",
      license             = "GPL",
      url                 = "http://desy.cfel.de/cid/cmi/cmistark",
      version             = "1.1-devel",
      long_description    = long_description,
      packages            = ['cmistark'],
      scripts             = ['scripts/cmistark_brute-force-orientation',
                             'scripts/cmistark_calculate_energy',
                             'scripts/cmistark_plot_energy',
                             'scripts/cmistark_print_energy'],
      test_suite          = 'tests',
      )
