.. _installing:

Installing pyXSIM
=================

Dependencies
------------

pyXSIM is compatible with Python 2.7 or 3.5+, and requires the following Python packages:

- `yt <http://yt-project.org>`_ (version 3.3.1 or higher)
- `NumPy <http://www.numpy.org>`_
- `AstroPy <http://www.astropy.org>`_
- `h5py <http://www.h5py.org>`_
- `six <https://pythonhosted.org/six/>`_

pyXSIM also has the following optional dependencies:

- `mpi4py <http://pythonhosted.org/mpi4py/>`_ (for running in parallel)

Installing
----------

pyXSIM can be installed in a few different ways. The simplest way is via the conda package if
you have the `Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_:

.. code-block:: bash

    [~]$ conda install -c jzuhone pyxsim

The second way to install pyXSIM is via pip. pip will attempt to download the dependencies and 
install them, if they are not already installed in your Python distribution:

.. code-block:: bash

    [~]$ pip install pyxsim

Alternatively, to install into your Python distribution from `source <http://github.com/jzuhone/pyxsim>`_:

.. code-block:: bash

    [~]$ python setup.py install
