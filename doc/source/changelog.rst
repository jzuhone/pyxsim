.. _changelog:

ChangeLog
=========

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