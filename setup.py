

#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

setup(name = 'pluto-crater-counting',
      description = 'Package for interactively counting Pluto and Charon craters',
      author = 'C.M. Gosmeyer',
      url = 'https://github.com/cgosmeyer/pluto-crater-counting',
      packages = find_packages(),
      install_requires = ['astropy', 'matplotlib', 'numpy', 'skimage']
     )

