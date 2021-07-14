.. _installing:

Installing pyXSIM
=================

Dependencies
------------

pyXSIM 3.x is compatible with Python 3.7+, and requires the following Python
packages:

- `yt <http://yt-project.org>`_ (version 4.0 or higher)
- `soxs <http://hea-www.cfa.harvard.edu/~jzuhone/soxs>`_ (version 3.0.0 or 
  higher)
- `NumPy <http://www.numpy.org>`_
- `SciPy <http://www.scipy.org>`_
- `AstroPy <http://www.astropy.org>`_
- `h5py <http://www.h5py.org>`_

pyXSIM also has the following optional dependencies:

- `mpi4py <http://pythonhosted.org/mpi4py/>`_ (for running in parallel)

Installing
----------

pyXSIM can be installed in a few different ways. The simplest way is via the 
conda package if you have the 
`Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_:

.. code-block:: bash

    [~]$ conda install -c jzuhone -c astropy pyxsim

Note both the ``jzuhone`` and ``astropy`` channels are required. This will install all of
the necessary dependencies.

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
