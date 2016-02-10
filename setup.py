#!/usr/bin/env python

from setuptools import setup, find_packages

entry_points = {
    'console_scripts': [
        'rat = rat.cli:dispatch',
    ]}



setup(name='rat',
      version='0.1',
      description='RNAseq Tools',
      entry_points=entry_points,
      include_package_data=True,
      packages=find_packages(),
      author='mark Fiers',
      author_email='mark.fiers.42@gmail.com',
      url=''
     )
