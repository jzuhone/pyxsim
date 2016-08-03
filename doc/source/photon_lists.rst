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
method to create a new :class:`~pyxsim.photon_list.PhotonList`. The arguments taken by 
:meth:`~pyxsim.photon_list.PhotonList.from_data_source` are as follows:

* ``data_source``: A :class:`~yt.data_objects.data_containers.YTSelectionContainer`. The 
  3D data source from which the photons will be generated.
* ``redshift``: The cosmological redshift for the photons. Determines the angular diameter
  and luminosity distances to the source given a ``cosmology``, which determines the number
  of photons. 
* ``area``: The collecting area to determine the number of photons. If units are
  not specified, it is assumed to be in :math:`cm^2`.
* ``exp_time``: The exposure time to determine the number of photons. If units are
  not specified, it is assumed to be in seconds.
* ``source_model`` : A :class:`~pyxsim.source_models.SourceModel` used to generate the
  photons. see :ref:`source-models` for options.
* ``parameters`` (optional): A dictionary of parameters to be passed to the source model if 
  necessary.
* ``center`` (optional): string or array-like object. The origin of the photon spatial 
  coordinates. Accepts "c", "max", or a coordinate. If not specified, pyXSIM attempts to use 
  the "center" field parameter of the data_source. 
* ``dist`` (optional): The angular diameter distance, only should be used for nearby sources. 
  This may be supplied instead of it being determined from the ``redshift`` and ``cosmology``.
  If units are not specified, it is assumed to be in Mpc. To use this, the redshift must be 
  set to zero. 
* ``cosmology`` (optional): A :class:`~yt.utilities.cosmology.Cosmology` object which supplies 
  cosmological informaton. If the ``data_source`` has cosmological parameters, they will be
  used, otherwise a :math:`\Lambda{\rm CDM}` cosmology with the following parameters are assumed: 
  :math:`H_0` = 71 km/s/Mpc, :math:`\Omega_m` = 0.27, :math:`\Omega_\Lambda` = 0.73. 
* ``velocity_fields`` (optional): A list of fields to use for the velocity to Doppler-shift the 
  photons. If not specified, the following will be assumed:   
  ``['velocity_x', 'velocity_y', 'velocity_z']`` for grid datasets, and 
  ``['particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z']`` for particle datasets.

As an example, we'll assume we have created a ``source_model`` representing the thermal emission 
from the plasma (see :ref:`source-models` for more details on how to create one): 

.. code-block:: python

    redshift = 0.05 # The redshift to the object. 
    area = (3000., "cm**2") # A constant effective area to generate the photons with.
    exp_time = (100., "ks") # The exposure time to generate the photons with. 
    center = sp.center # A center in 3D for the photon positions. If not specified, 
                       # the center of the `data_source` will be chosen.
    
    # Optionally, construct a cosmology object. 
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

:class:`~pyxsim.photon_list.PhotonList` instances which have been written to files can be
merged together, using the :func:`~pyxsim.utils.merge_files` function. This may be useful 
if you have so many photons to generate that they do not fit into memory all in one go.

:func:`~pyxsim.utils.merge_files` takes a list of input filenames, and an output filename. 
The optional keyword arguments are ``clobber``, which decides whether or not an existing file 
will be overwritten, and ``add_exposure_times`` decides whether or not the final file will 
have an exposure time of the sum of the times in the separate files or that of the longest 
exposure time between the files. 

.. code-block:: python

    from pyxsim import merge_files
    merge_files(["photons_0.h5","photons_1.h5","photons_3.h5"], "photons.h5",
                clobber=True, add_exposure_times=True)
