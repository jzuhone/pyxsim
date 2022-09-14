.. _source-models:

Source Models for X-ray Emission
================================

pyXSIM comes with a number of pre-defined ``SourceModel`` types for 
generating X-ray emission. Source models form the core of pyXSIM. 
These ``SourceModel``s can be used in three different ways:

* Generating X-ray emission fields for yt
* Generating spectra from yt data containers such as spheres, rectangular
  solids, etc.
* Generating mock photon lists for use in instrument simulators

There are currently three different types of source models in pyXSIM:

.. toctree::
    :maxdepth: 1

    thermal_sources
    powerlaw_sources
    line_sources
