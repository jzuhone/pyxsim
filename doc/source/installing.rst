Installing ``photon_simulator``
-------------------------------

Dependencies
++++++++++++

``photon_simulator`` is compatible with Python 2.7 or 3.5, and requires the following Python packages:

- `yt <http://yt-project.org>`_
- `NumPy <http://www.numpy.org>`_
- `AstroPy <http://www.astropy.org>`_
- `h5py <http://www.h5py.org>`_
- `six <https://pythonhosted.org/six/>`_

``photon_simulator`` also has the following optional dependencies:

- `mpi4py <http://pythonhosted.org/mpi4py/>`_ (for running in parallel)
- `PyXspec <http://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/>`_ (for generating certain types of spectral models)

Installing
++++++++++

``photon_simulator`` can be installed using pip. pip will attempt to download the dependencies and 
install them, if they are not already installed in your Python distribution. For an easy
installation of the dependencies, using a Python package distribution is recommended. For
example, using the `Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_:
  
.. code-block:: bash

    [~]$ conda install yt astropy h5py
    
Installing these three should be sufficient to get the remaining dependencies. In either case, install
``photon_simulator`` using pip:

.. code-block:: bash

    [~]$ pip install photon_simulator

Or, to install into your Python distribution from `source <http://bitbucket.org/jzuhone/photon_simulator>`_:

.. code-block:: bash

    [~]$ python setup.py install

Or, to install to a local directory, use:

.. code-block:: bash

    [~]$ python setup.py install --prefix=/path/to/location/

Then make sure your ``PYTHONPATH`` points to this location.

Using PyXspec for Spectral Models
+++++++++++++++++++++++++++++++++

For this to work, XSPEC and PyXspec must be compiled with the same Python distribution that 
``photon_simulator`` and its dependencies are in. Currently, PyXspec only works with Python 2.7. 

Issues on Mac OS X
==================

There may be additional issues with compiling and using ``photon_simulator`` and PyXspec on 
your Mac OS X system. 
