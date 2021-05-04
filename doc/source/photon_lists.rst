.. _photon-lists:

Photon Lists
============

The first step in creating synthetic X-ray observations with pyXSIM is that of
generating a photon list. A photon list is a three-dimensional HDF5 dataset that
contains a distribution of cell positions, velocities, and photon energies
within each cell. Once created, it be used to generated simulated X-ray events
from a particular line of sight. 

Generating a New Photon List from a Data Source
-----------------------------------------------

A photon list is generated from a data source, by which we mean any spatially
extended object in three dimensions represented in yt. For example, we may have
a dataset of a galaxy cluster simulation, and wish to generate photons from the
thermal emission of the hot cluster plasma. We can load up the dataset in yt and
create a spherical data source centered on the cluster potential minimum:

.. code-block:: python
    
    import yt
    from yt.utilities.cosmology import Cosmology
    import pyxsim
    
    ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")
    sp = ds.sphere("c", (1.0, "Mpc"))
    
This sphere object ``sp`` can then be used as input to the 
:func:`~pyxsim.photon_list.make_photons` function to create a new photon list.
The arguments taken by :func:`~pyxsim.photon_list.make_photons` are as follows:

* ``photon_prefix``: The prefix of the filename(s) to be written. If run in 
  serial, the filename will be ``"{photon_prefix}.h5"``, if run in parallel, the 
  filenames will be ``"{photon_prefix}.{mpi_rank}.h5"``.
* ``data_source``: A :class:`~yt.data_objects.data_containers.YTSelectionContainer`. 
  The 3D data source from which the photons will be generated.
* ``redshift``: The cosmological redshift for the photons. Determines the 
  angular diameter and luminosity distances to the source given a ``cosmology``,
  which determines the number of photons. 
* ``area``: The collecting area to determine the number of photons. If units are
  not specified, it is assumed to be in :math:`cm^2`.
* ``exp_time``: The exposure time to determine the number of photons. If units
  are not specified, it is assumed to be in seconds.
* ``source_model`` : A :class:`~pyxsim.source_models.SourceModel` used to 
  generate the photons. see :ref:`source-models` for options.
* ``parameters`` (optional): A dictionary of parameters to be passed to the 
  source model if necessary.
* ``point_sources`` (optional): If True, the photons will be assumed to be
  generated from the exact positions of the cells or particles and not smeared
  around within a volume. Default: False
* ``center`` (optional): string or array-like object. The origin of the photon
  spatial coordinates. Accepts "c", "max", or a coordinate. If not specified, 
  pyXSIM attempts to use the "center" field parameter of the data_source. 
* ``dist`` (optional): The angular diameter distance, only should be used for
  nearby sources. This may be supplied instead of it being determined from the 
  ``redshift`` and ``cosmology``. If units are not specified, it is assumed to
  be in Mpc. To use this, the redshift must be set to zero. 
* ``cosmology`` (optional): A :class:`~yt.utilities.cosmology.Cosmology` object
  which supplies cosmological informaton. If the ``data_source`` has 
  cosmological parameters, they will be used, otherwise a 
  :math:`\Lambda{\rm CDM}` cosmology with the following parameters are assumed: 
  :math:`H_0` = 71 km/s/Mpc, :math:`\Omega_m` = 0.27, 
  :math:`\Omega_\Lambda` = 0.73. 
* ``velocity_fields`` (optional): A list of fields to use for the velocity to
  Doppler-shift the photons. If not specified, the following will be assumed:   
  ``['velocity_x', 'velocity_y', 'velocity_z']`` for grid datasets, and 
  ``['particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z']`` 
  for particle datasets.

As an example, we'll assume we have created a ``source_model`` representing the
thermal emission from the plasma (see :ref:`source-models` for more details on
how to create one): 

.. code-block:: python

    redshift = 0.05 # The redshift to the object. 
    area = (3000., "cm**2") # A constant effective area to generate the photons with.
    exp_time = (100., "ks") # The exposure time to generate the photons with. 
    center = sp.center # A center in 3D for the photon positions. If not specified, 
                       # the center of the `data_source` will be chosen.
    
    # Optionally, construct a cosmology object. 
    cosmo = Cosmology(hubble_constant=0.68, omega_matter=0.31, omega_lambda=0.69)
    
    n_photons, n_cells = pyxsim.make_photons("my_photons", sp, redshift, area,
                                             exp_time, source_model, 
                                             center=center, cosmology=cosmo)

If you run on one core, this will write a file called ``"my_photons.h5"`` 
containing the photon list. If run on (say) 6 cores, it will write 6 files,
called ``"my_photons.[0-5].h5"``. The total number of photons is returned in
``n_photons``, and the total number of cells with photons is returned in
``n_cells``.

If you want to simulate photons from a a nearby object, set the redshift to zero
and specify a distance using the ``dist`` keyword argument:

.. code-block:: python

    n_photons, n_cells = pyxsim.make_photons("my_photons", sp, 0.0, area, 
                                             exp_time, source_model, 
                                             center=center, dist=(4., "kpc"))

By default, the photons generated from the cells or particles in the simulation 
will be smeared throughout the volume of those elements. To treat all of the 
cells or particles in the dataset as if they are point sources, set 
``point_sources=True``:

.. code-block:: python

    n_photons, n_cells = pyxsim.make_photons("my_photons", sp, redshift, area,
                                             exp_time, source_model, 
                                             center=center, point_sources=True)

By default, for computing the Doppler shifts of the photons, pyXSIM uses the 
default velocity fields of the dataset, which are ``"velocity_x"``, 
``"velocity_y"``, and ``"velocity_z"`` for grid/cell-based datasets and 
``"particle_velocity_x"``, ``"particle_velocity_y"``, and 
``"particle_velocity_z"`` for particle-based datasets. If you need to use other
fields, you can specify them using the ``velocity_fields`` keyword argument:

.. code-block:: python

    vfields = ["velx", "vely", "velz"]
    n_photons, n_cells = pyxsim.make_photons("my_photons", sp, redshift, area,
                                             exp_time, source_model, 
                                             center=center, dist=(4., "kpc"), 
                                             velocity_fields=vfields)


Merging Photon Lists
--------------------

Photon lists which have been written to files can be merged together, using the 
:func:`~pyxsim.utils.merge_files` function. This may be useful if you generate photons from
different sources or source types that are co-spatial.

:func:`~pyxsim.utils.merge_files` takes a list of input filenames, and an output filename. 
The optional keyword arguments are ``overwrite``, which decides whether or not an existing file 
will be overwritten, and ``add_exposure_times`` decides whether or not the final file will 
have an exposure time of the sum of the times in the separate files or that of the longest 
exposure time between the files. 

.. code-block:: python

    from pyxsim import merge_files
    merge_files(["photons_0.h5","photons_1.h5","photons_3.h5"], "photons.h5",
                overwrite=True, add_exposure_times=True)