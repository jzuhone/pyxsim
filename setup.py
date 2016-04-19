#!/usr/bin/env python
from setuptools import setup
from setuptools.extension import Extension
import numpy as np

cython_extensions = [
    Extension("pyxsim.cutils",
              sources=["pyxsim/cutils.pyx"],
              language="c", libraries=["m"],
              include_dirs=[np.get_include()])]

setup(name='pyxsim',
      packages=['pyxsim'],
      version='0.1.0',
      description='Python tools for ACIS Ops',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://bitbucket.org/jzuhone/pyxsim',
      install_requires=["six","numpy","astropy","h5py","yt"],
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      ext_modules=cython_extensions,
      )
