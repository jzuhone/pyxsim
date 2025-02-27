.. _changelog:

ChangeLog
=========

Version 4.4.2
-------------

This version of pyXSIM contains a bugfix for simulations where particle-based
emission sources have no intrinsic widths (neither cell sizes or smoothing
lengths).

Version 4.4.1
-------------

This version of pyXSIM contains a bugfix, where supplying an array of floats
for abundance tables in thermal source models was not working correctly.

Version 4.4.0
-------------

This version of pyXSIM contains a bugfix and updates compatibility with other
packages.

* Support for Python 3.8 has been dropped.
* The minimum required version of yt has been bumped up to 4.3.0.
* Cython 3.0 is now supported as a compile-time dependency for those building
  pyXSIM from source.
* A subtle bug that caused photons created from cells or particles on the edges
  of input regions to be scattered to large distances has been fixed.

Version 4.3.0
-------------

This version of pyXSIM contains new features.

* In conjunction with SOXS 4.6.0, it is now possible to pass an HDF5 events file
  created by :func:`~pyxsim.photon_list.project_photons` directly to the SOXS
  instrument simulator for use in event simulation. See :ref:`instr-soxs` for more
  details.
* The method :meth:`~pyxsim.sources.SourceModel.make_source_fields` now returns
  a fourth field, which is the X-ray "count rate" in photons/s. See :ref:`xray-fields`
  for more details.
* It is now possible to change the field names supplied by
  :meth:`~pyxsim.sources.SourceModel.make_source_fields` and
  :meth:`~pyxsim.sources.SourceModel.make_intensity_fields` if one desires to by
  supplying a value for the optional ``band_name`` keyword argument. If specified,
  this argument will replace the ``"{emin}_{emax}_keV`` part of the field name. See
  :ref:`xray-fields` for more details.
* Two new convenience methods, :meth:`~pyxsim.sources.SourceModel.make_line_source_fields`
  and :meth:`~pyxsim.sources.SourceModel.make_line_intensity_fields`, have been
  added, which create new fields for the line emission from a source. See :ref:`line-fields`
  for more details.
* For thermal emission sources, is is now possible to supply a ``min_entropy``
  parameter, which will set a minimum entropy for the gas cells or particles to
  be considered in the emission calculation. See :ref:`thermal-sources` for more
  details.

Version 4.2.1
-------------

This version of pyXSIM contains two bugfixes.

* A bug that sometimes prevented X-ray emission fields from being created
  in yt using thermal source models if a spatially varying metallicity field
  was used has been fixed.
* A bug that prevented luminosity field names from being displayed correctly in
  yt plots has been fixed.

Version 4.2.0
-------------

This version of pyXSIM contains one new feature and a critical bugfix.

* A bug introduced upstream in SOXS 4.4.0 caused normalization of emission and
  spectra in the :class:`~pyxsim.source_models.thermal_sources.IGMSourceModel`
  to be incorrect. This is fixed in SOXS 4.5.1 and this version requires that
  version of SOXS to work correctly.
* The function :func:`~pyxsim.photon_list.make_photons` now has a ``fields_to_keep``
  option that allows one to store fields other than the positions, velocities, and
  sizes of cells or particles to the photon list HDF5 file. See :ref:`generate_new`
  for more details.

Version 4.1.1
-------------

This version of pyXSIM fixes a bug in the
:class:`~pyxsim.source_models.thermal_sources.IGMSourceModel` where low density,
low temperature, and low metallicity plasmas would sometimes result in spectra with
negative values.

Version 4.1.0
-------------

This version of pyXSIM adds new features and improvements, as well as a bug fix.

* Installation and use on Windows 64-bit platforms is now supported.
* Generating photons from thermal emission models should be somewhat faster
  (sometimes up to 50% faster) thanks to more efficient spectral interpolation
  routines.
* A new class for reading information about and creating spectra from photon lists,
  :class:`~pyxsim.photon_list.PhotonList`, has been added. See :ref:`photon-list-class`
  for more details.
* Much more information about the parameters used to create photon lists and event
  lists are stored in the HDF5 files, as well as the versions of yt, pyXSIM, and SOXS
  used to create them. This info can be obtained from the
  :class:`~pyxsim.event_list.EventList` class, as well as the new
  :class:`~pyxsim.photon_list.PhotonList` class.
* The calculation of the Doppler shifting of photons from particle or cell velocities
  was incorrect for velocities in the relativistic regime and has been fixed. This did
  not affect anyone working in the Newtonian regime.
* All-sky projections can now utilize a ``center_vector`` which decides what position
  in the dataset the center of the sky coordinate system points to. See :ref:`allsky`
  for more details.
* The minimum version of yt has been bumped up to 4.1.3.
* The minimum version of SOXS has been bumped up to 4.3.0.
* This version of pyXSIM uses the new versions of the spectral models for the Cloudy CIE
  and IGM models provided in SOXS. See the information about them `here <http://hea-www.cfa.harvard.edu/soxs/users_guide/thermal_spectra.html#cloudy-cie-spectrum-generators>`_
  and `here <http://hea-www.cfa.harvard.edu/soxs/users_guide/thermal_spectra.html#igm-spectrum-generators>`_.



Version 4.0.1
-------------

This is a critical bugfix update to pyXSIM. The only change is that this version
has a minimum requirement of yt 4.1.2, which fixes a bug in which the number density
of hydrogen nuclei for certain Arepo datasets (including those from Illustris TNG)
was a factor of 2 too high. This, in turn, resulted in emission measures that were
a factor 2 too high, and also affects calculations that include photoionization,
which depend on the number density of hydrogen nuclei.

The bug, and its fix, are detailed here:

https://github.com/yt-project/yt/pull/4211

Version 4.0.0
-------------

This is a major update to pyXSIM which incldues many new features.

* In addition to the already existing mode of generating photon lists for
  producing synthetic observations, pyXSIM now has two new modes: creating
  X-ray emission, luminosity, and intensity fields for use in yt (see :ref:`xray-fields`)
  and creating spectra from yt data containers (see :ref:`xray-spectra`).
* Many changes to thermal sources, all of which are detailed in :ref:`thermal-sources`:
    * Modeling thermal emission spectra has been refactored into three new classes:
      :class:`~pyxsim.source_models.thermal_sources.CIESourceModel`,
      :class:`~pyxsim.source_models.thermal_sources.NEISourceModel`, and
      :class:`~pyxsim.source_models.thermal_sources.IGMSourceModel`. For sources
      in CIE, it is now possible to use SPEX, MeKaL, and Cloudy CIE models for
      the spectra. The new IGM model includes photoionization and resonant
      scattering off of the CXB.
    * There is now no default value for the ``max_density`` parameter in
      ``ThermalSourceModel`` instances (previously it was ``max_density=5e-25``).
    * ``ThermalSourceModel`` subclasses can now use log-spaced energy binning
      for the spectral model.
    * It is now possible to specify a yt field that allows for the hydrogen
      fraction to spatially vary for thermal sources.
    * Abundance tables from Feldman (1996) and Cloudy 17.03 have been added as options
      for specifying solar abundances.
* When creating a photon list, it is now possible to add a ``bulk_velocity``
  parameter, which will change the frame of reference for the velocity fields
  used in Doppler shifting of photons. See :ref:`generating-photon-lists` for more
  details.
* A new (and experimental) mode for creating all-sky mock observations using
  "internal" observers has been added to pyXSIM. See :ref:`allsky` for more details.
* A new option for using a flat-field projection when producing an event list has
  been added. See :ref:`event-lists` for more details.
* A new option for saving the line-of-sight coordinates when producing an event list
  has been added. See :ref:`event-lists` for more details.
* When producing an event list, it is now possible to change the abundance table
  assumed if one is using the "TBabs" model for foreground Galactic absorption. See
  :ref:`event-lists` for more details.

Version 3.0.1
-------------

This is a bugfix release to pyXSIM. All users are strongly encouraged to upgrade.

* A bug in the :class:`~pyxsim.source_models.ThermalSourceModel` which resulted
  in crashes when encountering a chunk of one or zero gas particles has been fixed.
* When using variable elements in :class:`~pyxsim.source_models.ThermalSourceModel`,
  if the ``Zmet`` argument for the remaining elements was a field and was a mass
  fraction, its conversion to solar units was computed incorrectly. This has now
  been fixed.

Version 3.0.0
-------------

This major update to pyXSIM contains a number of updates, including some
backwards-incompatible changes to the API. To figure out how to transfer
your code to version 3.x, please read :ref:`v2_to_v3`.

* A brand-new paradigm for generating photon lists and event lists has been
  created. In the new scheme, one does not create ``PhotonList`` and ``EventList``
  classes, but instead two functions, :func:`~pyxsim.photon_list.make_photons`
  and :func:`~pyxsim.photon_list.project_photons` are used to create photon lists
  and event lists which are stored on disk as they are made. This allows for very
  large photon lists and event lists to be created without holding them all in
  memory at once. For guidance on how to use the new functions, see
  :ref:`photon-lists` and :ref:`event-lists`.
* Support for Python 2.7 has been dropped in this version.
* The minimum supported yt version is now 4.0.
* The ``Zmet`` keyword argument to
  :class:`~pyxsim.source_models.ThermalSourceModel` has been changed to a required
  argument.
* The default minimum temperature ``kT_min`` for the
  :class:`~pyxsim.source_models.ThermalSourceModel` has been changed from 0.008
  keV to 0.025 keV.
* The ``max_density`` keyword argument to
  :class:`~pyxsim.source_models.ThermalSourceModel` has had its default value
  set to :math:`5 \times 10^{-25} g cm^{-3}`; previously it was ``None``.
* The X-ray binaries source generator has been dropped from pyXSIM.
* The background and point-source source generators have been removed, as this
  functionality can now be used within SOXS.

Version 2.3.1
-------------

This version contains bug fixes.

* Bugs were fixed to ensure compatibility with both yt 3.x and yt 4.0 (beta).
* A bug was fixed that resulted in odd behavior of the progress bars when in
  a Jupyter notebook.

Version 2.3.0
-------------

This version contains bug fixes and minor enhancements.

* This version supports ``h5py`` 3.x, which deprecated the use of accessing
  HDF5 dataset data using the ``.value`` attribute.
* This version supports both the ``yt`` 3.x series and the ``yt`` 4.0 beta
  version.
* Previous versions of pyXSIM scaled thermal emission by :math:`n_en_{H+}`,
  where :math:`n_{H+}` is the number density of free protons. However, the
  correct scaling is :math:`n_en_{H}`, where :math:`n_{H}` is the number
  density of hydrogen. This has been fixed.
* A bug which occurred when variable individual elements were used in the
  :class:`~pyxsim.source_models.ThermalSourceModel` has been fixed.
* The progress bar now updates correctly for generating sky positions when
  creating a :class:`~pyxsim.event_list.EventList`.
* Some minor speedups have been achieved in the
  :class:`~pyxsim.source_models.ThermalSourceModel` class.

Version 2.2.0
-------------

This version contains feature enhancements (with some backwards-incompatible
changes) and optimizations.

* The 2.2.x series of pyXSIM will be the last to support Python 2.7.
* Support for non-equilibrium ionization plasma emission using AtomDB has been
  added to pyXSIM. see :ref:`nei` for more details.
* The default AtomDB/APEC version for pyXSIM is now v3.0.9.
* The ability to change the redshift, collecting area, exposure time, or
  distance of the source when creating a :class:`~pyxsim.event_list.EventList`
  from :meth:`~pyxsim.photon_list.PhotonList.project_photons` has been removed.
  This was a little-used feature that was potentially confusing to users, and
  is mostly unnecessary given that the photon number will be reduced when
  convolving with any instrument simulators. This change also made the code
  simpler and resulted in optimizations. The related keyword arguments to
  :meth:`~pyxsim.photon_list.PhotonList.project_photons` will still be accepted,
  but will be ignored.
* Arepo data is now fully supported.
* A new option to treat each cell or particle which emits photons as a point
  source has been added to the :meth:`~pyxsim.photon_list.PhotonList.from_data_source`
  method of :class:`~pyxsim.photon_list.PhotonList`.
* The built-in instrument models are now deprecated, as well as
  :class:`~pyxsim.event_list.ConvolvedEventList` objects. For convolution with
  instrument models, users are encouraged to use
  `SOXS <http://hea-www.cfa.harvard.edu/~jzuhone/soxs>`_ or another instrument
  simulation package.

Version 2.1.1
-------------

This version contains a single bugfix. The conversion factors between mass fractions and
solar units for individual elements in the :class:`~pyxsim.source_models.ThermalSourceModel`
were not being calculated correctly and has now been fixed. Simulations which used a single
metallicity field only were not affected by this bug.

Version 2.1.0
-------------

This version contains bugfixes and feature enhancements, as well new version requirements
for dependencies.

* This version of pyXSIM requires AstroPy version 2.0 or higher, yt version 3.4 or higher,
  and SOXS version 2.0 or higher.
* A number of bugs in the :func:`~pyxsim.utils.merge_files` function were fixed.
* The ``"redshift"`` and ``"d_a"`` parameters have been removed from
  :class:`~pyxsim.event_list.EventList` objects, as events at different redshifts/distances
  should be able to be combined together.
* If two :class:`~pyxsim.event_list.EventList` objects are added and their ``"sky_center"``
  parameters differ, the two :class:`~pyxsim.event_list.EventList` objects are added together and
  the ``"sky_center"`` parameter of the first one is used. Previously, two different
  ``"sky_center"`` parameters would have thrown an error.
* With the introduction of instrument models for ACIS-S in SOXS v2.0, it is no longer
  necessary to retain the ACIS-S response file with pyXSIM and in general response files
  will no longer be included with pyXSIM for instrument simulation.
* The ``ACIS_I`` and ``ACIS_S`` instrument models have been updated from Cycle 18 to Cycle 19.
* The ability to use separate abundances of individual elements in the computation of
  a thermal spectrum has been added to the :class:`~pyxsim.source_models.ThermalSourceModel`.
  See :ref:`thermal-sources` and :ref:`var-abund` for more information.
* In the creation of a :class:`~pyxsim.source_models.ThermalSourceModel`, it is now possible
  to use Solar abundance tables other than the implicitly assumed Anders & Grevesse 1989. See
  and :ref:`thermal-sources` and :ref:`solar-abund-tables` for details.
* It is now possible to simulate a :class:`~pyxsim.source_models.ThermalSourceModel` without
  emission lines. See :ref:`thermal-sources` for details.
* :meth:`~pyxsim.photon_list.PhotonList.project_photons` has been refactored under the hood
  to improve memory usage and speed.

Version 2.0.0
-------------

This is a major new release of pyXSIM, which fixes bugs, adds a number of new features,
but most importantly, implements a simpler API in many aspects. A number of the changes
in this version are backwards-incompatible with previous versions, and where applicable
is noted below. A useful summary of the API changes with some code examples can be
found at :ref:`v1_to_v2`.

The largest (and largely hidden) change in this release is the outsourcing of
much of pyXSIM's capabilities to `SOXS <http://hea-www.cfa.harvard.edu/~jzuhone/soxs>`_,
which is a spin-off package from pyXSIM which models thermal spectra, foreground
galactic absorption, and convolving with instrument models. This results in far
less duplication between the code bases of these two closely related projects.

New features:

* A new class, :class:`~pyxsim.light_cone.XrayLightCone`, has been added which takes
  a number of redshift snapshots from a cosmological simulation and produces a light
  cone simulation of events from them. This is an experimental feature which should
  be considered in "beta", and currently only works with Enzo or Gadget-based
  cosmological simulations.
* A module has been added to generate X-ray photons from a population of X-ray
  binaries, both low-mass and high-mass. This assumes as input a simulation with star
  particles which have masses, ages, and metallicities. See :ref:`xray-binaries` for
  more information. This is an experimental feature which should be considered in "beta".
* A minor feature, but methods and functions that accept arguments such as ``area`` and
  ``exp_time`` which accept values with unit information can now accept
  :class:`~astropy.units.Quantity` instances.

Changes related to thermal source modeling:

* pyXSIM now uses SOXS to implement APEC-based thermal spectral models.
* The previously deprecated XSPEC-based thermal spectral models have been
  completely removed from this version, as they proved too difficult to maintain.
* It is no longer necessary to create a thermal spectral model object explicitly,
  as this is now handled by :class:`~pyxsim.source_models.ThermalSourceModel`.
  This method now takes the name of the spectral model as a parameter. Consequently,
  arguments needed for the creation of spectra now need to be passed to
  :class:`~pyxsim.source_models.ThermalSourceModel` upon creation of a new instance.
  This is a backwards-incompatible change.
* Thermal broadening of spectral lines is now on by default.

Changes related to modeling of foreground Galactic absorption:

* pyXSIM now uses SOXS to implement the `wabs` and `tbabs` foreground absorption
  models.
* The previously deprecated XSPEC-based spectral absorption models have been
  completely removed from this version, as they proved too difficult to maintain.
* It is no longer necessary to create a spectral absorption model object explicitly,
  as this is now handled by :meth:`~pyxsim.photon_list.PhotonList.project_photons`.
  This method now takes the name of the absorption model as a parameter. Consequently,
  the ``nH`` parameter for the hydrogen column is now a parameter which is passed
  to :meth:`~pyxsim.photon_list.PhotonList.project_photons`. This is a
  backwards-incompatible change.

The following changes arise from a refactor of ``InstrumentSimulator``

* The ``InstrumentSimulator`` class now uses the SOXS machinery for convolving with
  instrumental responses.
* The only operations performed by ``InstrumentSimulator`` are convolution with the
  effective area curve (using the ARF) and with the response matrix (using the RMF).
  No spatial PSF convolutions or rebinning operations can be applied. For more detailed
  instrument simulation, users are advised to write events to SIMPUT files and use SOXS directly.
* New *Hitomi* response files have been supplied with this version.
* The ``XRS_Imager`` and ``XRS_Calorimeter`` instruments have been renamed to
  ``Lynx_Imager`` and ``Lynx_Calorimeter``.

The following interrelated changes arise from a refactor of :class:`~pyxsim.event_list.EventList`:

* Instrument simulators now return a new :class:`~pyxsim.event_list.ConvolvedEventList`
  instance, which contains the data and parameters for convolved events. It is no longer
  possible for :class:`~pyxsim.event_list.EventList` instances to contain convolved events.
* The :meth:`~pyxsim.event_list.EventList.write_spectrum` now only bins on unconvolved
  energy (see next bullet for the new way to bin on channel).
* The new :class:`~pyxsim.event_list.ConvolvedEventList` class has a method,
  :meth:`~pyxsim.event_list.ConvolvedEventList.write_channel_spectrum`, which writes a
  spectrum binned on PI or PHA channels.
* :class:`~pyxsim.event_list.EventList` instances no longer contain pixelated coordinates
  for events based on the resolution of the simulation, but only sky coordinates. The
  :meth:`~pyxsim.event_list.EventList.write_fits_file` and
  :meth:`~pyxsim.event_list.EventList.write_fits_image` methods now accept arguments
  which create custom pixelizations for event files and images.
* :class:`~pyxsim.event_list.EventList` instances no longer contain all events on all
  processors when created in parallel, but each processor now contains a subset of the
  events. The I/O routines for :class:`~pyxsim.event_list.EventList` have been rewritten
  so that all events are still written to the file.
* The methods for generating events from point sources and backgrounds have been removed
  from :class:`~pyxsim.event_list.EventList` and now exist as "source generators" which
  return new event lists. See :ref:`source-generators` for more information.

Other changes:

* The ``sky_center`` parameter to :meth:`~pyxsim.photon_list.PhotonList.project_photons`
  is now a required argument. This is a backwards-incompatible change.
* The ``clobber`` keyword argument for overwriting files has been changed to ``overwrite``.
  This is a backwards-incompatible change.
* Handling for `cut regions <http://yt-project.org/doc/analyzing/filtering.html#cut-regions>`_
  when creating a :class:`~pyxsim.photon_list.PhotonList` for a dataset with periodic
  boundaries has been improved in this release.
* :class:`~pyxsim.photon_list.PhotonList` and :class:`~pyxsim.event_list.EventList`
  instances now use the same keys as their corresponding HDF5 files. The old keys will
  still work for the time being, but are deprecated. This is a backwards-incompatible
  change.
* The optional argument ``smooth_positions`` has been added to the
  :meth:`~pyxsim.photon_list.PhotonList.project_photons` method, which allows one to
  smooth the event positions to avoid block-shaped artifcats in images with lots of
  counts.
* Thermal spectral models no longer require a ``cleanup_spectrum`` method. Spectral
  absorption models no longer require ``setup_spectrum`` and ``cleanup_spectrum``
  methods. Source models no longer require a ``cleanup_model`` method.
* pyXSIM now has `SciPy <http://www.scipy.org>`_ as a required dependence.
* Throughout the code, pseudo-random number generators can now be specified simply
  as integer seeds in signatures to functions which take the keyword argument ``prng``.

Version 1.2.6
-------------

This is a bugfix release that ensures that fields with units of ``code_metallicity`` are
properly handled.

Version 1.2.5
-------------

This is a bugfix release with two fixes:

* Ensured that metallicity fields in the :class:`~pyxsim.source_models.ThermalSourceModel`
  are properly scaled to the Anders & Grevasse (1989) solar metallicity since this is
  what APEC assumes.
* Support for octree mesh datasets (such as RAMSES) has now been added.

Version 1.2.4
-------------

This version fixes a single bug, ensuring that the metallicity is converted to
solar units in thermal source models.

Version 1.2.3
-------------

This is a bugfix release.

* Gadget binary (non-HDF5) datasets are now supported.
* Make sure that SPH datasets assume fully ionized gas if an ``ElectronAbundance`` field is not present.
* The normalization of the power-law and line emission models was incorrect by a factor of :math:`1/(1+z)`.
  This has been fixed.

Version 1.2.2
-------------

This is a bugfix release.

* Position fields for SPH datasets will now be correctly detected for
  irregularly shaped sources.
* Photon numbers for all sources are now being generated assuming a Poisson
  distribution.
* pyXSIM will no longer automatically emit a deprecation warning when it tries
  to import ``assert_same_wcs`` from yt.
* Minor documentation fixes.

Version 1.2.1
-------------

This is a bugfix release.

* Fixed a bug when writing FITS table files when AstroPy 1.3 is installed.
* Fixed an import error which occurs when using the yt development branch.
* Minor documentation updates

Version 1.2.0
-------------

This version contains bugfixes and performance enhancements, as well as a new test suite.

* We are now running a test suite which automatically checks changes to the code pushed up to the
  `GitHub repository <http://github.com/jzuhone/pyxsim>`_.
* The definition of the ``norm`` parameter for the :meth:`~pyxsim.spectral_models.TableApecModel.return_spectrum`
  method is now consistent with the `normal Xspec definition <http://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelApec.html>`_.
* Annoying NumPy indexing warnings have been silenced by only using signed ints for indexing.
* Absorption models have been refactored to have a more common structure.
* For table-based absorption models, the cross-section is now interpolated instead of the absorption factor itself,
  which should be more accurate.
* XSpec-based spectral models are officially in deprecation; they will be removed in a future release.
* A bug that prevented response matrices from not being read properly with old versions of AstroPy was fixed.

Version 1.1.1
-------------

This version is a bugfix and optimization release.

* Some speedups have been achieved in the convolution of energies with RMFs.
* An error is now thrown if one attempts to use a zero or negative redshift in
  :meth:`~pyxsim.photon_list.PhotonList.from_data_source` without specifying a distance.

Version 1.1.0
-------------

This version contains a bugfix and some minor new features.

* Fixed a bug which did not use the correct file names for AtomDB tables when using
  ``TableApecModel``.
* Refactored the absorption model handling into a new class. No user-facing changes have been made.
* Added special classes for the TBabs and wabs absorption models.
* De-emphasizing XSpec-based spectral models in favor of the table-based alternatives.

Version 1.0.1
-------------

This is solely a bugfix release.

* Ensured that spherical and box-shaped regions which wrap periodic boundaries are
  handled correctly.
* The width of event list field of view is determined correctly for 3-D source
  distributions with high aspect ratios.
