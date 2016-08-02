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
  when the :class:`~pyxsim.photon_list.PhotonList` was created.
* ``exp_time_new`` (optional): The exposure time to assume for the observation. Used to reduce the number
  of events from the initially large sample of photons. The default value is the value used when the 
  :class:`~pyxsim.photon_list.PhotonList` was created.

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

If convolved with responses using an instrument model (see :ref:`instruments` for more
details), these FITS files are "standard" events files which may be read and analyzed 
by other X-ray software tools such as ds9, CIAO, HEATOOLS, etc.

SIMPUT
++++++

An :class:`~pyxsim.event_list.EventList` can be exported to the SIMPUT file format for
reading in by other packages that simulate particular instruments, such as MARX or SIMX
(see :ref:`instruments` for more details). This is done by calling the 
:meth:`~pyxsim.event_list.EventList.write_simput_file` method:

.. code-block:: python

    events.write_simput_file("my_great_events", clobber=False, emin=0.1, emax=9.0)

where the first argument is the prefix for the files that will be created (the SIMPUT 
file and a photon list sidecar file), and the other optional arguments control whether
or not an existing file will be overwritten and the minimum and maximum energies of the
events written to the file. Currently, SIMPUT files are used for export only; they
cannot be used to read events back into pyXSIM. 

Adding Background and Point Source Events
-----------------------------------------

Methods are provided for adding background and point source events to an existing 
:class:`~pyxsim.event_list.EventList`. To add background photons, call
:meth:`~pyxsim.event_list.EventList.add_background`:

.. code-block:: python

    from numpy.random import RandomState
    prng = RandomState(25)

    events.add_background(ebins, spec, prng=prng, absorb_model=tbabs_model)


Manipulating Event Lists
------------------------


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

    events.write_spectrum()
