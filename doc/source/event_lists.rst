.. _event-lists:

Event Lists
===========

Creating a New Event List by Projecting from a Photon List
----------------------------------------------------------

Adding Background and Point Source Events
-----------------------------------------

Manipulating Event Lists
------------------------

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

These FITS files are "standard" events files which may be read by other X-ray 
software tools such as ds9, CIAO, etc.

SIMPUT
++++++

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


