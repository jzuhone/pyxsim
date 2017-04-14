.. _event-lists:

Event Lists
===========

:class:`~pyxsim.event_list.EventList` objects are the main data products of pyXSIM, since
they represent the synthetic observations that may be compared to and analyzed in the same
fashion as real observed data. :class:`~pyxsim.event_list.EventList`\s can be produced from
:class:`~pyxsim.photon_list.PhotonList`\s by projecting along a particular axis, or read in
from disk from a previous projection. They can also be manipulated in a number of ways,
combined with other sets of events, and used to produce other useful data products. 

Creating a New Event List by Projecting from a Photon List
----------------------------------------------------------

New :class:`~pyxsim.event_list.EventList`\s are created from a :class:`~pyxsim.photon_list.PhotonList`
using the :meth:`~pyxsim.photon_list.PhotonList.project_photons` method. This method accepts a
line of sight direction as its only required argument, with a number of additional keyword 
arguments to customize the resulting projection. The arguments are:

* ``normal``: The line of sight direction to project along. Accepts either a coordinate axis (``"x"``,
  ``"y"``, or ``"z"``), or a three-vector for an off-axis projection, e.g. ``[1.0, -0.3, 0.24]``. 
* ``area_new`` (optional): The (constant) collecting area to assume for the observation. Used to reduce
  the number of events from the initially large sample of photons. The default value is the value used 
  when the :class:`~pyxsim.photon_list.PhotonList` was created. Units are in :math:`cm^2`.
* ``exp_time_new`` (optional): The exposure time to assume for the observation. Used to reduce the number
  of events from the initially large sample of photons. The default value is the value used when the 
  :class:`~pyxsim.photon_list.PhotonList` was created. Units are in seconds.
* ``redshift_new`` (optional): The value of the redshift to assume for the observation. Used to reduce the
  of events from the initially large sample of photons. The default value is the value used when the 
  :class:`~pyxsim.photon_list.PhotonList` was created.
* ``dist_new`` (optional): The value of the angular diameter to assume for the observation. Use for nearby
  sources instead of the redshift. If units are not specified, it is assumed to be in Mpc. Used to reduce the
  of events from the initially large sample of photons. The default value is the value used when the 
  :class:`~pyxsim.photon_list.PhotonList` was created. To use this, the redshift must be set to zero. 
* ``absorb_model`` (optional): :class:`~pyxsim.spectral_models.AbsorptionModel` A model for foreground 
  galactic absorption. See :ref:`absorb-models` for details on how to construct one.  
* ``sky_center`` (optional): Central RA, Dec of the events in degrees. Default ``(30.0, 45.0)``.
* ``no_shifting`` (optional): If set to True, the photon energies will not be velocity Doppler shifted. Default False.
* ``north_vector`` (optional): A vector defining the "up" direction, e.g. ``[0.0, 1.0, 0.0]``.
  This option sets the orientation of the plane of projection. If not set, an arbitrary grid-aligned 
  ``north_vector`` is chosen. Ignored in the case where a particular axis (e.g., "x", "y", or "z") is 
  explicitly specified.
* ``prng`` (optional): A pseudo-random number generator, :class:`~numpy.random.RandomState` object, or
  :mod:`~numpy.random` is the default. Use this if you have a reason to generate the same set of random 
  numbers, such as for a test. 

This is also the stage where foreground galactic absorption can be applied. See :ref:`absorb-models` for
details on how to construct models for absorption. 

Assuming one then has a :class:`~pyxsim.photon_list.PhotonList` ``photons``, example invocations could look
like this:

A simple projection along an axis:

.. code-block:: python

    events = photons.project_photons("z")
        
An off-axis projection with altered exposure time and redshift:

.. code-block:: python

    events = photons.project_photons([0.1, -0.3, 0.5], area_new=(200., "cm**2"), redshift_new=1.0)

An on-axis projection with absorption and given a particular coordinate center:

.. code-block:: python

    abs_model = pyxsim.TBabsModel(0.01)
    events = photons.project_photons("y", absorb_model=abs_model, sky_center=(12.0, -30.0))

An off-axis projection with a ``north_vector``, without Doppler velocity shifting, and a specific
random number generator:

.. code-block:: python
    
    from numpy.random import RandomState
    prng = RandomState(25)
    events = photons.project_photons([0.1, -0.3, 0.5], no_shifting=True, north_vector=[1.0,0.0,0.0],
                                     prng=prng)

.. note::

    Unlike the ``photon_simulator`` analysis module in yt, the ability to convolve the event energies
    using an ARF and RMF has been taken out of this step entirely and moved into a new instrument 
    simulator step. See :ref:`instruments` for details. 
    
Saving/Reading Raw Events to/from Disk
--------------------------------------

For storage and later usage, events can be written to disk and read back in later
in three file formats. 

HDF5
++++

Any :class:`~pyxsim.event_list.EventList` instance may be saved to disk in the
convenient HDF5 file format by calling the :meth:`~pyxsim.event_list.EventList.write_h5_file`
method:

.. code-block:: python
    
    events.write_h5_file("cluster_events.h5")
    
To read previously stored events back from disk, use the 
:meth:`~pyxsim.event_list.EventList.from_h5_file` method:

.. code-block:: python

    events = EventList.from_h5_file("cluster_events.h5")

FITS
++++

Any :class:`~pyxsim.event_list.EventList` instance may be saved to disk in the
FITS format by calling the :meth:`~pyxsim.event_list.EventList.write_fits_file`
method:

.. code-block:: python

    events.write_fits_file("cluster_events.fits", clobber=True)
    
The ``clobber`` keyword argument is used to allow (or prevent) overwrites of 
files if they already exist. To read previously stored events back from disk, 
use the :meth:`~pyxsim.event_list.EventList.from_fits_file` method:

.. code-block:: python

    events = EventList.from_fits_file("cluster_events.fits")

.. note::

    If the events have convolved with responses using a built-in instrument model (see 
    :ref:`instruments` for more details), then the characteristics of the convolved events
    are stored in the HDF5 or FITS file, including the spectral channels. In the case of the
    FITS files, they are "standard" events files which may be read and analyzed by other X-ray 
    software tools such as ds9, CIAO, HEATOOLS, etc.

.. _simput:

SIMPUT
++++++

An :class:`~pyxsim.event_list.EventList` can be exported to the SIMPUT file format for
reading in by other packages that simulate particular instruments, such as
`MARX <http://space.mit.edu/ASC/MARX/>`_ or `SIMX <http://hea-www.cfa.harvard.edu/simx/>`_
(see also :ref:`instruments`). This is done by calling the 
:meth:`~pyxsim.event_list.EventList.write_simput_file` method:

.. code-block:: python

    events.write_simput_file("my_great_events", clobber=False, emin=0.1, emax=9.0)

where the first argument is the prefix for the files that will be created (the SIMPUT 
file and a photon list sidecar file), and the other optional arguments control whether
or not an existing file will be overwritten and the minimum and maximum energies of the
events written to the file. Currently, SIMPUT files are used for export only; they
cannot be used to read events back into pyXSIM. 

Manipulating Event Lists
------------------------

There are a couple of options for manipulating :class:`~pyxsim.event_list.EventList` objects. 

If two :class:`~pyxsim.event_list.EventList` objects were created with the same parameters (e.g.
exposure time, collecting area, etc.), and only the events are different, they can be simply added
together to return a new :class:`~pyxsim.event_list.EventList`:

.. code-block:: python

    events = events1 + events2
    
An error will be thrown if the parameters do not match between the two lists. 

The second way an :class:`~pyxsim.event_list.EventList` can be changed is by using a region file.
This requires the `pyregion <http://pyregion.readthedocs.io/>`_ package to be installed. If you have
a region file, simply provide it to the :meth:`~pyxsim.event_list.EventList.filter_events` method:

.. code-block:: python

    some_events = events.filter_events("annulus.reg")

which creates a new :class:`~pyxsim.event_list.EventList` object with only the events which fall within
the region. 

Saving Derived Products from Event Lists
----------------------------------------

:class:`~pyxsim.event_list.EventList` instances can produce binned images and spectra
from their events. Both products are written in FITS format.

Images
++++++

To produce a binned image, call the :meth:`~pyxsim.event_list.EventList.write_fits_image`
method:

.. code-block:: python

    events.write_fits_image("myimage.fits", clobber=True, emin=0.5, emax=7.0)

which writes an image binned at the finest resolution of the simulation in the file
``"myimage.fits"``. Set ``clobber=True`` if the file is already there and you 
want to overwrite it. The ``emin`` and ``emax`` parameters control the energy range
of the events which will be included in the image (default is to include all of the
events).

Spectra
+++++++

To produce a binned spectrum, call :meth:`~pyxsim.event_list.EventList.write_spectrum`. 

.. code-block:: python

    events.write_spectrum("myspec.fits", bin_type="energy", emin=0.1,
                          emax=10.0, nchan=2000, clobber=False)

In this case, we have chosen ``"energy"`` as the ``bin_type``, which will bin the unconvolved
event energies using the ``emin``, ``emax``, and ``nchan`` arguments into a histogram which
will be written to the file as a spectrum. If ``bin_type="channel``:

.. code-block:: python

    events.write_spectrum("myspec.fits", bin_type="channel")

the spectrum will be binned by channel instead. This is particularly useful if the events
were convolved with an instrument simulator, as this spectrum can be read and fit using programs
such as XSPEC. As usual, the ``clobber`` argument determines whether or not a file can
be overwritten. 