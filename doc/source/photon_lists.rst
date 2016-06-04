.. _photon-lists:

Photon Lists
============

The first step in creating synthetic X-ray observations with pyXSIM is that of generating
a :class:`~pyxsim.photon_list.PhotonList`. The :class:`~pyxsim.photon_list.PhotonList` is 
a three-dimensional object that represents a distribution of cell positions, velocities,
and photon energies within each cell. Once created or restored from disk, a 
:class:`~pyxsim.photon_list.PhotonList` can be used to generated simulated X-ray events
from a particular line of sight. 

Generating Photons from a Data Source
-------------------------------------

The typical starting point for a :class:`~pyxsim.photon_list.PhotonList` is to generate
one from a data source, by which we mean any spatially extended object in three dimensions
represented in yt. To demonstrate how to create a :class:`~pyxsim.photon_list.PhotonList`,
we'll use a 

Saving/Reading Photons to/from Disk
-----------------------------------

Any :class:`~pyxsim.photon_list.PhotonList` instance may be saved to disk in the convenient
HDF5 format by calling the :meth:`~pyxsim.photon_list.PhotonList.write_h5_file` method:

.. code-block:: python
    
    photons.write_h5_file("cluster_photons.h5")
    
This writes the photon positions, velocities, length scales, energies, and associated
parameters to disk. To read previously stored photons back from disk, use the
:meth:`~pyxsim.photon_list.PhotonList.from_file` method:

.. code-block:: python

    photons = PhotonList.from_file("cluster_photons.h5")

Merging Photon Lists
--------------------


