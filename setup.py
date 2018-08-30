#!/usr/bin/env python
from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy as np
import versioneer

cython_extensions = [
    Extension("pyxsim.lib.sky_functions",
              ["pyxsim/lib/sky_functions.pyx"],
              language="c", libraries=["m"],
              include_dirs=[np.get_include()])
]

setup(name='pyxsim',
      packages=find_packages(),
      version=versioneer.get_version(),
      description='Python package for simulating X-ray observations of astrophysical sources',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/jzuhone/pyxsim',
      setup_requires=["numpy", "cython>=0.24"],
      install_requires=["six","numpy","astropy>=2.0","h5py","scipy","yt>=3.4","soxs>=2.1"],
      include_package_data=True,
      ext_modules=cython_extensions,
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      )
