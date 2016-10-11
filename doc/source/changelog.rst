.. _changelog:

ChangeLog
=========

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