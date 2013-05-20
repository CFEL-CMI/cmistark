#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008,2009,2012 Jochen Küpper <software@jochen-kuepper.de>


import os
from distutils.core import setup

extra_compile_args = []
library_dirs = []

long_description = """JK Python extensions for Stark effect calculations

Python extensions for Stark effect calculations of molecules.

Original authors:   Jochen Küpper <software@jochen-kuepper.de>
Current maintainer: Jochen Küpper <software@jochen-kuepper.de>
"""


setup(name="cmistark",
      author              = "Jochen Küpper",
      author_email        = "software@jochen-kuepper.de",
      description         = "JK Python Stark effect extensions",
      license             = "GPL",
      url                 = "http://desy.cfel.de/cid/cmi",
      version             = "0.1.0",
      long_description    = long_description,
      package_dir         = {'cmistark': 'lib'},
      packages            = ['cmistark'],
      scripts             = ['scripts/cmistark_brute-force-orientation',
                             'scripts/cmistark_calculate_energy',
                             'scripts/cmistark_plot_energy',
                             'scripts/cmistark_print_energy']
      )
