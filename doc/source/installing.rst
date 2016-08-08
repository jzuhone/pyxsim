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
- `PyXspec <http://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/>`_ (for generating certain types of spectral models)

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

Using PyXspec for Spectral Models
---------------------------------

pyXSIM provides the capability of using thermal and foreground absorption models from
`XSPEC <https://heasarc.gsfc.nasa.gov/xanadu/xspec/>`_ using the
`PyXspec <https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/>`_ package. For this to
work, you must compile XSPEC and PyXspec with the same Python distribution that
pyXSIM and its dependencies are in. This can be done by simply compiling in an environment where
the correct ``python`` executable is the first in the path. For more details on compiling XSPEC,
visit the `"Installing HEASoft" <http://heasarc.gsfc.nasa.gov/docs/software/lheasoft/install.html>`_ page.
Currently, PyXspec only works with Python 2.7.

There are some issues you may encounter with generating spectral models using PyXspec. They are:

Conflicts with AstroPy
++++++++++++++++++++++

Both AstroPy and PyXspec use the CFITSIO library. Unfortunately, both of them also compile it as 
a shared library themselves and link to it. This can result in conflicts between the two. You may see
error messages like this:

.. code-block:: console

    ERROR: Mismatch in the CFITSIO_SONAME value in the fitsio.h include file
    that was used to build the CFITSIO library, and the value in the include file
    that was used when compiling the application program:
        Version used to build the CFITSIO library   = 1
        Version included by the application program = 5

    Fix this by recompiling and then relinking this application program
    with the CFITSIO library.
    
The only way to resolve this problem at present is to avoid importing the ``xspec`` and ``astropy``
modules in the same Python environment, script, or notebook. 