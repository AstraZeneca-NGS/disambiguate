#!/usr/bin/env python
from distutils.core import setup

version = "1.0"

setup(name='disambiguate',
      version=version,
      description='Script to disambiguate reads mapping to multiple genomes.',
      author='Miika Ahdesmaki',
      license="MIT",
      scripts=['disambiguate.py'])
