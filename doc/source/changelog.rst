.. _changelog:

ChangeLog
=========

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