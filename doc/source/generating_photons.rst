Generating Photon Lists
=======================

The first step in creating synthetic X-ray observations with ``photon_simulator`` is 
that of generating a ``PhotonList``. The 

Generating Photons from a Data Source
-------------------------------------

Saving/Reading Photons to/from Disk
-----------------------------------

Any ``PhotonList`` instance may be saved to disk by calling the ``write_h5_file`` method:

.. code-block:: python
    
    photons.write_h5_file("cluster_photons.h5")
    
This writes the photon positions, velocities, length scales, energies, and associated
parameters to disk in the convenient HDF5 file format. To read previously stored photons
back from disk, use the ``from_file`` method of the ``PhotonList`` class:

.. code-block:: python

    photons = PhotonList.from_file("cluster_photons.h5")

Merging Photon Lists
--------------------


