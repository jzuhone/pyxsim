.. _how-it-works:

How pyXSIM Works
================

This section details the physical models and basic algorithms of pyXSIM. 

Assumptions
-----------


From Source to Detector: Three Basic Steps
------------------------------------------

Step 1: Generate Photons
++++++++++++++++++++++++

We begin with a three-dimensional model of a source, whether from a dataset from some type of dynamical
simulation or from a dataset constructed by hand. This dataset is comprised of cells or particles, and 
is understandable by yt. Examples of this may be:
 
* A model of a galaxy cluster plasma, where the emission arises from various interactions between 
  the thermal electrons and ions
* A model of relativistic particles which emit inverse-Compton emission in a power-law spectrum 
* A model of dark matter particles where line emission arises from annihilation or decay of the particles

In any case, there is a fair amount of flexibility in the type of source and its emission that can be
modeled. 

Step 2: Project Photons
+++++++++++++++++++++++

Once we have generated photons in three dimensions, we can create a set of simulated events by projecting 
them along a particular line of sight. 

Step 3: Convolve with Instrumental Responses
++++++++++++++++++++++++++++++++++++++++++++

The implementation of this step is described in :ref:`instruments`. 

Limitations
-----------

Given the above framework, pyXSIM currently has the following limitations, which may or may not be 
lifted in the future: 

* pyXSIM is currently unable to simulate any sources with any optical thickness, absorption, or scattering.
* pyXSIM is currently unable to simulate any sources with explicit time dependence.
