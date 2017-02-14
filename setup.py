#!/usr/bin/env python
from setuptools import setup

version = "1.0"

setup(name='disambiguate',
      version=version,
      description='Script to disambiguate reads mapping to multiple genomes.',
      author='Miika Ahdesmaki',
      license="MIT",
      install_requires = [
          'pysam>=0.8.4',
          ],
      py_modules = [
          'disambiguate',
          ],
      entry_points = {
          'console_scripts': [
              'disambiguate = disambiguate:main',
              ],
          },
      )
