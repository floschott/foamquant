#!/usr/bin/env python

from distutils.core import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(name='foamquant',
      version='0.1.0',
      description='3D quantification of cellular material structure and dynamics',
      author='Scientific Data, MAX IV Laboratory, Lund University',
      author_email='florian_schott@outlook.fr',
      url='https://www.lth.se',
      packages=['foamquant'],
      package_dir={'foamquant': 'docs/source/FoamQuant'},
      install_requires=required,
     )
