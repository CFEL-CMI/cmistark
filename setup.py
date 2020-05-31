#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# Copyright (C) 2008,2009,2012,2013,2014,2015,2018,2020 Jochen Küpper <jochen.kuepper@cfel.de>


import os
from setuptools import setup

extra_compile_args = []
library_dirs = []

copyright = '2008,2014,2020 Jochen Küpper and the CFEL Controlledc Molecule Imaging group'
name = "cmistark"
version = "1.4.dev0"
release = version
long_description = """CMI Python extensions for Stark effect calculations

Python extensions for Stark effect calculations of molecules.
Developed by Jochen Küpper and the Controlled Molecule Imaging Group at the Center for Free-Electron Laser Science,
Hamburg, Germany.

Original author:    Jochen Küpper <jochen.kuepper@fhi-berlin.mpg.de>
Current maintainer: Jochen Küpper <jochen.kuepper@cfel.de>
"""


setup(name="cmistark",
      python_requires     = '>=3.5',
      author              = "Jochen Küpper and the CFEL Controlled Molecule Imaging group",
      author_email        = "jochen.kuepper@cfel.de",
      maintainer          = "Jochen Küpper",
      maintainer_email    = "jochen.kuepper@cfel.de",
      url                 = "https://www.controlled-molecule-imaging.org/research/further_projects/software/cmistark",
      description         = "CMI Python Stark effect calculations",
      version             = version,
      long_description    = long_description,
      license             = "GPLv3",
      packages            = ['cmistark'],
      scripts             = ['scripts/cmistark_brute-force-orientation',
                             'scripts/cmistark_calculate_energy',
                             'scripts/cmistark_plot_energy',
                             'scripts/cmistark_print_energy'],
      command_options={
          'build_sphinx': {
              'project': ('setup.py', name),
              'version': ('setup.py', version),
              'release': ('setup.py', release),
              'source_dir': ('setup.py', 'doc'),
              'copyright': ('setup.py', copyright)}
      },
      test_suite          = 'tests',
      classifiers=[
          'Development Status :: 6 - Mature',
          'Environment :: Console',
          'Environment :: MacOS X',
          'Intended Audience :: Education',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Science/Research',
          'License :: DFSG approved',
          'License :: OSI Approved',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX',
          'Operating System :: POSIX :: Linux',
          'Operating System :: Unix',
          'Operating System :: Other OS',
          'Programming Language :: Python :: 3',
          'Topic :: Education',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Physics'
          ],
      )
