.. _photon-lists:

Photon Lists
============

The first step in creating synthetic X-ray observations with pyXSIM is that of generating
a :class:`~pyxsim.photon_list.PhotonList`. The :class:`~pyxsim.photon_list.PhotonList` is 
a three-dimensional object that represents a distribution of cell positions, velocities,
and photon energies within each cell. Once created or restored from disk, a 
:class:`~pyxsim.photon_list.PhotonList` can be used to generated simulated X-ray events
from a particular line of sight. 

Generating a New Photon List from a Data Source
-----------------------------------------------

The typical starting point for a :class:`~pyxsim.photon_list.PhotonList` is to generate
one from a data source, by which we mean any spatially extended object in three dimensions
represented in yt. For example, we may have a dataset of a galaxy cluster simulation, and wish
to generate photons from the thermal emission of the hot cluster plasma. We can load up the
dataset in yt and create a spherical data source centered on the cluster potential minimum:

.. code-block:: python
    
    import yt
    from yt.utilities.cosmology import Cosmology
    import pyxsim
    
    ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")
    sp = ds.sphere("c", (1.0, "Mpc"))
    
This sphere object ``sp`` can then be used as input to the :meth:`~pyxsim.photon_list.PhotonList.from_data_source`
method to create a new :class:`~pyxsim.photon_list.PhotonList`. 
:meth:`~pyxsim.photon_list.PhotonList.from_data_source` requires a ``source_model`` argument; you
can find out how to setup a source model in the section on :ref:`source-models`. For now, we'll
assume we have created the ``source_model`` representing the thermal emission from the plasma. 

.. code-block:: python

    redshift = 0.05 # The redshift to the object. 
    area = (3000., "cm**2") # A constant effective area to generate the photons with.
    exp_time = (100., "ks") # The exposure time to generate the photons with. 
    center = sp.center # A center in 3D for the photon positions. Can be anything. If
                       # not specified, the center of the `data_source` will be chosen.
    
    # Optionally, construct a cosmology object. If you don't, then one will be chosen
    # with hubble_constant=0.71, omega_matter=0.27, omega_lambda=0.73.
    cosmo = Cosmology(hubble_constant=0.68, omega_matter=0.31, omega_lambda=0.69)
    
    photons = pyxsim.PhotonList.from_data_source(sp, redshift, area, exp_time,
                                                 source_model, center=center, 
                                                 cosmology=cosmo)

If you want to simulate photons from a a nearby object, set the redshift to zero
and specify a distance using the ``dist`` keyword argument:

.. code-block:: python

    photons = pyxsim.PhotonList.from_data_source(sp, 0.0, area, exp_time,
                                                 source_model, center=center, 
                                                 dist=(4., "kpc"))

By default, for computing the Doppler shifts of the photons, pyXSIM uses the default velocity 
fields of the dataset, which are ``"velocity_x"``, ``"velocity_y"``, and ``"velocity_z"`` 
for grid-based datasets and ``"particle_velocity_x"``, ``"particle_velocity_y"``, and 
``"particle_velocity_z"`` for particle-based datasets. If you need to use other fields, you 
can specify them using the ``velocity_fields`` keyword argument:

.. code-block:: python

    photons = pyxsim.PhotonList.from_data_source(sp, 0.0, area, exp_time,
                                                 source_model, center=center, 
                                                 dist=(4., "kpc"), 
                                                 velocity_fields=["velx", "vely", "velz"])

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


