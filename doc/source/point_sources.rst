.. _point-sources:

Point Sources
=============

The :func:`~pyxsim.source_generators.point_sources.make_point_sources` function is provided
for creating an :class:`~pyxsim.event_list.EventList` composed of events from point sources.

:func:`~pyxsim.source_generators.point_sources.make_point_sources` takes the following arguments:

* ``area``: float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`.
  The collecting area to determine the number of events. If units are not
  specified, it is assumed to be in :math:`\rm{cm}^2`.
* ``exp_time``: float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`.
  The exposure time to determine the number of events. If units are not specified, 
  it is assumed to be in seconds.
* ``positions``: array of source positions, shape 2xN. The positions of the point 
  sources in RA, Dec, where N is the number of point sources. Coordinates should 
  be in degrees.
* ``sky_center``: array-like. Center RA, Dec of the events in degrees.
* ``spectra``: list (size N) of :class:`~soxs.spectra.Spectrum` objects. The 
  spectra for the point sources, where N is the number of point sources. 
  Assumed to be in the observer frame.
* ``prng``: integer or :class:`~numpy.random.RandomState` object, optional.
  A pseudo-random number generator. Typically will only be specified if you 
  have a reason to generate the same set of random numbers, such as for a test. 
  Default is to use the :mod:`numpy.random` module.

A simple example: 

.. code-block:: python

    from yt import YTArray, YTQuantity
    import numpy as np
    from numpy.random import RandomState
    from soxs import Spectrum

    # Simulate two point sources with power-law spectra

    alpha1 = -1.0 # The power-law slopes of the two point sources' spectra
    alpha2 = -1.5
    
    redshift1 = 0.01 # The redshifts of the two point sources
    redshift2 = 0.02
    
    norm1 = 1.0e-10 # The normalizations of the two point sources' spectra 
    norm2 = 1.0e-9  # units are photons/s/cm**2/keV
        
    spec1 = Spectrum.from_powerlaw(alpha1, redshift1, norm1,
                                   emin=0.01, emax=20.0, nbins=10000)
    spec2 = Spectrum.from_powerlaw(alpha2, redshift2, norm2,
                                   emin=0.01, emax=20.0, nbins=10000)
    
    # Apply foreground absorption to the spectra
    
    spec1.apply_foreground_absorption(0.02)
    spec2.apply_foreground_absorption(0.02)

    # The RA, Dec positions of the two sources in degrees
    
    positions = [(30.01, 44.99), (29.98, 45.03)]
    
    area = (20000.0, "cm**2")
    exp_time = (400.0, "ks")
    sky_center = (10.0, 12.0) # in degrees

    events = make_point_sources(area, exp_time, positions, sky_center, 
                                [spec1, spec2], prng=23)

.. note::

    If you want to have absorbed spectra for the point sources, this should 
    be done on the :class:`~soxs.spectra.Spectrum` objects first using the 
    :meth:`~soxs.spectra.Spectrum.apply_foreground_absorption` method.
