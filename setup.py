#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name='pyxsim',
      packages=find_packages(),
      version='2.0.0',
      description='Python package for simulating X-ray observations of astrophysical sources',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/jzuhone/pyxsim',
      setup_requires=["numpy"],
      install_requires=["six","numpy","astropy>=1.3","h5py","scipy","yt>=3.3.5","soxs>=1.2.0"],
      include_package_data=True,
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      )
