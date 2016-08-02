.. _how-it-works:

How pyXSIM Works
================

This section details in brief the model, assumptions, and basic algorithms of pyXSIM.

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
modeled. The basic idea is that for each cell in a specified region of the source, photons are generated
assuming the physical properties of the call and an emission model for X-rays based on those physical 
properties. The number of photons to be generated is controlled by specifying a fiducial exposure time, 
collecting area, and redshift and/or distance to the source. Typically, at this stage the values of these
parameters are chosen such that they generate a large number of photons, much larger than would be actually
observed, since the goal is to create an initial distribution that will be a statistically robust sample
from which to draw events when creating synthetic observations. The photons are generated with three-dimensional
positions in the rest frame of the source (except that their energies are potentially cosmologically redshifted
with the given redshift), since the positions on the sky and the Doppler shifting of the photons are 
taken into account in the next step when a line of sight is provided. 

Step 2: Project Photons
+++++++++++++++++++++++

Once we have generated our large sample of photons in the previous step, we can create a set of simulated 
events by projecting them along a particular line of sight. 

Step 3: Convolve with Instrumental Responses
++++++++++++++++++++++++++++++++++++++++++++

Finally, the last step is to take the sky-projected event positions and energies and convolve them with
instrumental responses to produce realistic images and spectra that can be processed with standard software
tools for working with X-ray observations. pyXSIM provides a very simplified approach to convolving with
instrumental responses, but also provides a way to export the simulated events for use with other packages
which simulate X-ray observatories more accurately.

The details on how to convolve with instrumental responses can be found in :ref:`instruments`. 

Limitations
-----------

Given the above framework, pyXSIM currently has the following limitations, which may or may not be 
lifted in the future: 

* pyXSIM is currently unable to simulate any sources with any optical thickness, absorption, or scattering.
* pyXSIM is currently unable to simulate any sources with explicit time dependence.
* pyXSIM is currently limited to datasets in Cartesian coordinates. 
