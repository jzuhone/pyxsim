.. pyxsim documentation master file, created by
   sphinx-quickstart on Sat Apr  9 10:00:48 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyXSIM Documentation
====================

What is pyXSIM?
---------------

The Heritage of pyXSIM
----------------------

pyXSIM is an implementation of the `PHOX <http://www.mpa-garching.mpg.de/~kdolag/Phox/>`_
algorithm, developed for constructing mock X-ray observations from SPH datasets by
Veronica Biffi and Klaus Dolag. There are two relevant papers:

`Biffi, V., Dolag, K., Bohringer, H., & Lemson, G. 2012, MNRAS, 420,
3545 <http://adsabs.harvard.edu/abs/2012MNRAS.420.3545B>`_

`Biffi, V., Dolag, K., Bohringer, H. 2013, MNRAS, 428,
1395 <http://adsabs.harvard.edu/abs/2013MNRAS.428.1395B>`_

pyXSIM had a previous life as the ``photon_simulator`` analysis module as a part of the
`yt Project <http://yt-project.org>`_. pyXSIM still depends critically on yt to provide the
link between the simulation data and the algorithm for generating the X-ray photons. For
detailed information about the design of the algorithm in yt, check out
`the SciPy 2014 Proceedings. <http://conference.scipy.org/proceedings/scipy2014/zuhone.html>`_

.. toctree::
   :maxdepth: 2

   physics
   installing
   generating_photons
   source_models
   generating_events


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

