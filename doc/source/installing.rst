.. _installing:

Installing pyXSIM
=================

Dependencies
------------

pyXSIM 4.x is compatible with Python 3.9+, and requires the following Python
packages:

- `yt <https://yt-project.org>`_ (version 4.1.3 or higher)
- `soxs <http://hea-www.cfa.harvard.edu/soxs>`_ (version 4.3.0 or
  higher)
- `NumPy <https://www.numpy.org>`_
- `SciPy <https://www.scipy.org>`_
- `AstroPy <https://www.astropy.org>`_
- `h5py <https://www.h5py.org>`_

pyXSIM also has the following optional dependencies:

- `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_ (for running in parallel)

Installing
----------

pyXSIM can be installed in a few different ways. The simplest way is via the
conda package if you have the
`Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_:

.. code-block:: bash

    [~]$ conda install -c conda-forge pyxsim

This will install all of the necessary dependencies.

The second way to install pyXSIM is via pip. pip will attempt to download the
dependencies and install them, if they are not already installed in your Python
distribution:

.. code-block:: bash

    [~]$ pip install pyxsim

Alternatively, to install into your Python distribution from
`source <http://github.com/jzuhone/pyxsim>`_:

.. code-block:: bash

    [~]$ git clone http://github.com/jzuhone/pyxsim
    [~]$ cd pyxsim
    [~]$ pip install .
