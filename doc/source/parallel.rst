.. _parallel:

Running in Parallel
===================

For large jobs, pyXSIM can be run in parallel using MPI. To do this, you need 
to have an MPI library, and the `mpi4py <http://mpi4py.readthedocs.io/>`_ 
package installed into your Python stack. To run scripts in parallel, import 
``yt`` before importing ``pyxsim``, and call ``yt.enable_parallelism()``
before doing anything else:

.. code-block:: python

    import yt
    yt.enable_parallelism()
    import pyxsim
    
    # rest of code goes here...

When pyXSIM is run in parallel, photon lists and event lists will be created
across multiple cores and split across multiple files when written to disk. 
