.. _generating-photon-lists:

Generating Photon Lists
=======================

A photon list in pyXSIM is a three-dimensional HDF5 dataset that contains a
distribution of cell or particle positions, velocities, and photon energies
within each cell or particle. Once created, it be used to generated simulated
X-ray events from a particular line of sight.

.. _generate_new:

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
* ``bulk_velocity`` (optional): A 3-element array or list specifying the local
  velocity frame of reference. If not a :class:`~yt.units.yt_array.YTArray`,
  it is assumed to have units of km/s. Default: [0.0, 0.0, 0.0] km/s.
* ``observer`` (optional): Whether the observer is "external", and the flat-sky
  approximation can be assumed for the purposes of projecting photons onto the
  sky, or the observer is "internal" and a spherical projection is assumed.
  Default is "external". For more details on the "internal" observing mode see
  :ref:`allsky`.
* ``fields_to_keep`` (optional): A list of fields to add to the HDF5 photon
  list in addition to the cell or particle positions, velocities, and sizes.
  Default: None

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

    vfields = [("flash", "velx"), ("flash", "vely"), ("flash", "velz")]
    n_photons, n_cells = pyxsim.make_photons("my_photons", sp, redshift, area,
                                             exp_time, source_model,
                                             center=center, dist=(4., "kpc"),
                                             velocity_fields=vfields)

We can also add other fields to the file using the ``fields_to_keep`` option,
which will not be used for photon projection in later steps but can be used for
diagnostic and/or analysis purposes. These represent the fields at the positions
where photons have been generated from. For example, to add density and
temperature fields to the file:

.. code-block:: python

    fields_to_keep = [("gas", "density"), ("gas", "temperature")]
    n_photons, n_cells = pyxsim.make_photons("my_photons", sp, redshift, area,
                                             exp_time, source_model,
                                             center=center, dist=(4., "kpc"),
                                             fields_to_keep=fields_to_keep)


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


.. _photon-list-class:

The ``PhotonList`` Class
------------------------

Generated photons stored to an HDF5 file can be examined with the
:class:`~pyxsim.photon_list.PhotonList` class. Currently, this class can be
used to examine parameters of a photon list in the ``parameters`` dictionary:

.. code-block:: python

    photons = pyxsim.PhotonList("therm_photons.h5")
    print(photons.parameters)

.. code-block:: pycon

    {'bulk_velocity': array([0., 0., 0.]),
     'center': array([0., 0., 0.]),
     'data_type': 'cells',
     'fid_area': 500.0,
     'fid_d_a': 122.21820987067642,
     'fid_exp_time': 100000.0,
     'fid_redshift': 0.03,
     'hubble': 0.71,
     'observer': 'external',
     'omega_lambda': 0.73,
     'omega_matter': 0.27,
     'velocity_fields': array([[b'gas', b'velocity_x'],
            [b'gas', b'velocity_y'],
            [b'gas', b'velocity_z']], dtype='|S10')}

as well as other pertinent information in the ``info`` dictionary:

.. code-block:: python

    print(photons.info)

.. code-block:: pycon

    {'data_source': 'YTSphere (UniformGridData): , center=[0. 0. 0.] cm, radius=1.5428387904811624e+24 cm',
     'dataset': 'UniformGridData',
     'pyxsim_version': '4.1b1.dev29+g1c09873.d20221228',
     'source_model': "CIESourceModel(
                          model=apec
                          emin=1 keV
                          emax=80.0 keV
                          nbins=5000
                          Zmet=0.3
                          binscale=linear
                          temperature_field=('gas', 'temperature')
                          emission_measure_field=('gas', 'emission_measure')
                          kT_min=0.025
                          kT_max=64.0
                          method=invert_cdf
                          model_vers=3.0.9
                          max_density=None
                          abund_table=angr
                          h_fraction=0.7065215023571868
                          var_elem={}
                          nolines=False
                          thermal_broad=True
                      )",
     'soxs_version': '4.2.1',
     'yt_version': '4.2.dev0'}

If this photon list file has originated from merged photon lists, then there
will be multiple instances of each piece of information, numbered by the
file, e.g. ``"soxs_version_0"``, ``"soxs_version_1"``, and so on. The original
files used to make the merge will be stored in the key ``"original_files"``.

Spectra
=======

To produce a spectrum binned on energy, call
:meth:`~pyxsim.photon_list.PhotonList.write_spectrum`.

.. code-block:: python

    specfile = "myspec.fits" # filename to write to
    emin = 0.1 # minimum energy of spectrum
    emax = 10.0 # maximum energy of spectrum
    nchan = 2000 # number of bins in spectrum
    photons.write_spectrum(specfile, emin, emax, nchan, overwrite=False)

This bins the photon energies using the ``emin``, ``emax``, and
``nchan`` arguments into a histogram which will be written to the file as a
spectrum. As usual, the ``overwrite`` argument determines whether or not a file
can be overwritten.
