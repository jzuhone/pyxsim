.. pyxsim documentation master file, created by
   sphinx-quickstart on Sat Apr  9 10:00:48 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyXSIM Documentation
====================

.. raw:: html

   <figure style="display: table; float: right; margin: 0 0 20px 20px;">
   <img src="_images/cluster_merger_events.png" width="400" style="float: right;"/>
   <figcaption style="display: table-caption; caption-side: bottom;">
   Simulated X-ray events from a binary cluster merger with mass ratio of 1:3, 
   initial impact parameter of 0.5 Mpc, at an epoch of 2.0 Gyr. Events were 
   convolved with the <em>Chandra</em> ACIS-I responses. Image taken from the 
   <a href="http://gcmc.hub.yt">Galaxy Cluster Merger Catalog</a>.
   </figcaption>
   </figure>

What is pyXSIM?
---------------

pyXSIM is a Python package for simulating X-ray observations from astrophysical 
sources.

X-rays probe the high-energy universe, from hot galaxy clusters to compact 
objects such as neutron stars and black holes and many interesting sources in 
between. pyXSIM makes it possible to generate synthetic X-ray observations of 
these sources from a wide variety of models, whether from grid-based simulation 
codes such as FLASH, Enzo, and Athena, to particle-based codes such as Gadget 
and AREPO, and even from datasets that have been created "by hand", such as from
NumPy arrays. pyXSIM also provides facilities for manipulating the synthetic 
observations it produces in various ways, as well as ways to export the 
simulated X-ray events to other software packages to simulate the end products
of specific X-ray observatories. 

The Heritage of pyXSIM
----------------------

pyXSIM is an implementation of the 
`PHOX <http://www.mpa-garching.mpg.de/~kdolag/Phox/>`_
algorithm, developed for constructing mock X-ray observations from SPH datasets
by Veronica Biffi and Klaus Dolag. There are two relevant papers:

`Biffi, V., Dolag, K., Bohringer, H., & Lemson, G. 2012, MNRAS, 420,
3545 <http://adsabs.harvard.edu/abs/2012MNRAS.420.3545B>`_

`Biffi, V., Dolag, K., Bohringer, H. 2013, MNRAS, 428,
1395 <http://adsabs.harvard.edu/abs/2013MNRAS.428.1395B>`_

pyXSIM had a previous life as the ``photon_simulator`` analysis module as a part
of the `yt Project <http://yt-project.org>`_. pyXSIM still depends critically on
yt to provide the link between the simulation data and the algorithm for 
generating the X-ray photons. For detailed information about the design of the 
algorithm in yt, check out
`the SciPy 2014 Proceedings. <http://conference.scipy.org/proceedings/scipy2014/zuhone.html>`_

License
-------

pyXSIM is released under a 
`BSD 3-clause license <https://opensource.org/licenses/BSD-3-Clause>`_.

Current Version
---------------

The current stable version is 3.0.0. See the :ref:`changelog` for details on 
changes from previous versions.

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   installing
   v2_to_v3
   how_it_works
   basic_yt_concepts
   photon_lists
   source_models
   absorb_models    
   event_lists
   light_cone
   instruments
   parallel
   getting_help
   cookbook/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

