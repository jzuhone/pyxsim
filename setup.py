#!/usr/bin/env python
import os

import numpy as np
from setuptools import find_packages, setup
from setuptools.extension import Extension

if os.name == "nt":
    std_libs = []
else:
    std_libs = ["m"]

cython_extensions = [
    Extension(
        "pyxsim.lib.sky_functions",
        ["pyxsim/lib/sky_functions.pyx"],
        language="c",
        libraries=std_libs,
        include_dirs=[np.get_include()],
    ),
    Extension(
        "pyxsim.lib.spectra",
        ["pyxsim/lib/spectra.pyx"],
        language="c",
        libraries=std_libs,
        include_dirs=[np.get_include()],
    ),
    Extension(
        "pyxsim.lib.interpolate",
        ["pyxsim/lib/interpolate.pyx"],
        language="c",
        libraries=std_libs,
        include_dirs=[np.get_include()],
    ),
]


setup(
    packages=find_packages(),
    url="http://github.com/jzuhone/pyxsim",
    include_package_data=True,
    ext_modules=cython_extensions,
)
