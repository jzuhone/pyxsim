.. _changelog:

ChangeLog
=========

Version 2.0.0
-------------

This is a major new release of pyXSIM, which fixes bugs, adds a number of new features,
but most importantly, implements a simpler API in many aspects. A number of the changes 
in this version are backwards-incompatible with previous versions, and where applicable
is noted below.

The largest (and largely hidden) change in this release is the outsourcing of 
much of pyXSIM's capabilities to [SOXS](http://hea-www.cfa.harvard.edu/~jzuhone/soxs), 
which is a spin-off package from pyXSIM which models thermal spectra, foreground
galactic absorption, and convolving with instrument models. This results in far 
less duplication between the code bases of these two closely related projects.

New features:

* A new class, :class:`~pyxsim.light_cone.XrayLightCone`, has been added which takes
  a number of redshift snapshots from a cosmological simulation and produces a light
  cone simulation of events from them. This is an experimental feature which should
  be considered in "beta".

Changes related to thermal source modeling:

* pyXSIM now uses SOXS to implement APEC-based thermal spectral models.
* The previously deprecated XSPEC-based thermal spectral models have been 
  completely removed from this version, as they proved too difficult to maintain. 
* It is no longer necessary to create a thermal spectral model object explicitly,
  as this is now handled by :class:`~pyxsim.source_models.ThermalSourceModel`.
  This method now takes the name of the spectral model as a parameter. Consequently, 
  arguments needed for the creation of spectra now need to be passed to 
  :class:`~pyxsim.source_models.ThermalSourceModel` upon creation of a new instance. 
* Thermal broadening of spectral lines is now on by default.

Changes related to modeling of foreground Galactic absorption:

* pyXSIM now uses SOXS to implement the `wabs` foreground absorption model.
* The previously deprecated XSPEC-based spectral absorption models have been 
  completely removed from this version, as they proved too difficult to maintain. 
* It is no longer necessary to create a spectral absorption model object explicitly,
  as this is now handled by :meth:`~pyxsim.photon_list.PhotonList.project_photons`.
  This method now takes the name of the absorption model as a parameter. Consequently, 
  the ``nH`` parameter for the hydrogen column is now a parameter which is passed 
  to :meth:`~pyxsim.photon_list.PhotonList.project_photons`.

The following changes arise from a refactor of :class:`~pyxsim.instruments.InstrumentSimulator`:

* The :class:`~pyxsim.instruments.InstrumentSimulator` class now uses the SOXS machinery
  for convolving with instrumental responses.
* The only operations performed by :class:`~pyxsim.instruments.InstrumentSimulator` are
  convolution with the effective area curve (using the ARF) and with the response matrix
  (using the RMF). No spatial PSF convolutions or rebinning operations can be applied. For
  more detailed instrument simulation, users are advised to write events to SIMPUT files
  and use SOXS directly. 
* New *Hitomi* response files have been supplied with this version. 

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

Other changes:

* The ``sky_center`` parameter to :meth:`~pyxsim.photon_list.PhotonList.project_photons`
  is now a required argument. This is a backwards-incompatible change.
* The ``clobber`` keyword argument for overwriting files has been changed to ``overwrite``.
  This is a backwards-incompatible change.
* :class:`~pyxsim.photon_list.PhotonList` and :class:`~pyxsim.event_list.EventList`
  instances now use the same keys as their corresponding HDF5 files. The old keys will 
  still work for the time being, but are deprecated. This is a backwards-incompatible 
  change.
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