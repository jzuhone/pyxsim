#!/usr/bin/env python
from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy as np

cython_extensions = [
    Extension("pyxsim.lib.sky_functions",
              ["pyxsim/lib/sky_functions.pyx"],
              language="c", libraries=["m"],
              include_dirs=[np.get_include()])
]


setup(name='pyxsim',
      packages=find_packages(),
      description='Python package for simulating X-ray observations of astrophysical sources',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/jzuhone/pyxsim',
      install_requires=["numpy", "astropy>=4.0", "h5py>=3.0", "scipy", "yt>=4.0.3", 
                        "soxs>=4.0b2", "tqdm"],
      include_package_data=True,
      ext_modules=cython_extensions,
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      )
