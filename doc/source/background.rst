.. _background:

Background
==========

The :func:`~pyxsim.source_generators.background.make_background` function is provided
for creating an :class:`~pyxsim.event_list.EventList` composed of events from an 
astrophysical background or foreground.

:func:`~pyxsim.source_generators.background.make_background` takes the following arguments:

* ``area``: float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`.
  The collecting area to determine the number of events. If units are not 
  specified, it is assumed to be in :math:`\rm{cm}^2`.
* ``exp_time``: float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`/
  The exposure time to determine the number of events. If units are not 
  specified, it is assumed to be in seconds.
* ``fov``: float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`.
  The field of view of the event file. If units are not provided, they are 
  assumed to be in arcminutes.
* ``sky_center``: array-like. Center RA, Dec of the events in degrees.
* ``spectrum``: :class:`~soxs.spectra.Spectrum` The spectrum for the background.
* ``prng``: integer or :class:`~numpy.random.RandomState` object, optional.
  A pseudo-random number generator. Typically will only be specified if you 
  have a reason to generate the same set of random numbers, such as for a test. 
  Default is to use the :mod:`numpy.random` module.

The background is assumed to be spatially constant over the entire field of view. 
A simple example: 

.. code-block:: python

    from soxs import Spectrum

    # Not terribly realistic, but this example simulates a constant background
    # with respect to energy. The input spectrum value here is in units of 
    # photons/s/cm**2/keV and will be spread out over the entire field of view.

    area = (20000.0, "cm**2")
    exp_time = (400.0, "ks")
    spec = Spectrum.from_constant(1.0e-5, emin=0.1, emax=10.0)
    sky_center = (21.0, 22.0) # in degrees
    fov = (20.0, "arcmin")
    
    events = make_background(area, exp_time, fov, sky_center, spec, prng=25)

.. note::

    If you want to have an absorbed spectrum for the background, this should 
    be done on the :class:`~soxs.spectra.Spectrum` object first using the 
    :meth:`~soxs.spectra.Spectrum.apply_foreground_absorption` method.

.. note::

    This method does not work for particle/instrumental backgrounds. That should
    be handled using your preferred specific instrument simulation package.