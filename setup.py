#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# Copyright (C) 2008,2009,2012,2013 Jochen Küpper <jochen.kuepper@cfel.de>


import os
from distutils.core import setup

extra_compile_args = []
library_dirs = []

long_description = """CMI Python extensions for Stark effect calculations

Python extensions for Stark effect calculations of molecules.
Developed by the Controlled Molecule Imaging Group at CFEL

Original author:     Jochen Küpper <jochen.kuepper@cfel.de>
Current maintainers: Yuan-Pin Chang <yuan.pin.chang@cfel.de> and Jochen Küpper <jochen.kuepper@cfel.de>
"""


setup(name="cmistark",
      author              = "Jochen Küpper, CFEL-CMI group, et al (see AUTHORS)",
      author_email        = "jochen.kuepper@cfel.de",
      description         = "CMI Python Stark effect extensions",
      license             = "GPL",
      url                 = "http://desy.cfel.de/cid/cmi",
      version             = "1.0",
      long_description    = long_description,
      package_dir         = {'cmistark': 'lib'},
      packages            = ['cmistark'],
      scripts             = ['scripts/cmistark_brute-force-orientation',
                             'scripts/cmistark_calculate_energy',
                             'scripts/cmistark_plot_energy',
                             'scripts/cmistark_print_energy']
      )
