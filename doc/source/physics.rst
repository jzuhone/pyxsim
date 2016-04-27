.. _physics:

The Physics Behind pyXSIM
=========================

If you don't care to know or already know how X-ray emission from astrophysical sources works, and just want to get
started, you can skip to :ref:`installing`.

Assumptions
-----------

From Source to Detector: Three Basic Steps
------------------------------------------

Step 1: Generate Photons
++++++++++++++++++++++++

Step 2: Project Photons
+++++++++++++++++++++++

Once we have generated photons in three dimensions, we can create a set of simulated events

Step 3: Convolve with Instrumental Responses
++++++++++++++++++++++++++++++++++++++++++++

Limitations
-----------

Given the above framework, pyXSIM currently has the following limitations, which may or may not be 
lifted in the future: 

* pyXSIM is currently unable to simulate any sources with any optical thickness, absorption, or scattering.
* pyXSIM is currently unable to simulate any sources with explicit time dependence.
