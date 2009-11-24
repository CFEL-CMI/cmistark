#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2008,2009 Jochen Küpper <software@jochen-kuepper.de>


import os
from distutils.core import setup

extra_compile_args = []
library_dirs = []

long_description = """JK Python extensions for Stark effect calculations

Python extensions for Stark effect calculations of molecules.

Original authors:   Jochen Küpper <software@jochen-kuepper.de>
Current maintainer: Jochen Küpper <software@jochen-kuepper.de>
"""


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
      scripts             = ['scripts/jkstark_brute-force-orientation',
                             'scripts/jkstark_calculate_energy',
                             'scripts/jkstark_plot_energy',
                             'scripts/jkstark_print_energy']
      )

