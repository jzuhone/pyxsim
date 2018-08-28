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
* ``sky_center``: Central RA, Dec of the events in degrees.
* ``absorb_model`` (optional): A string or :class:`~pyxsim.spectral_models.AbsorptionModel` class 
  representing a model for foreground galactic absorption. This parameter can take a string or the 
  class itself. See :ref:`absorb-models` for more details on how to use them. Known options for 
  strings are ``"wabs"`` and ``"tbabs"``.
* ``nH`` (optional): The foreground galactic column density in units of 
  :math:`10^{22} \rm{atoms} \rm{cm}^{-2}`, for use when one is applying foreground galactic absorption.
* ``no_shifting`` (optional): If set to True, the photon energies will not be velocity Doppler shifted. Default False.
* ``north_vector`` (optional): A vector defining the "up" direction, e.g. ``[0.0, 1.0, 0.0]``.
  This option sets the orientation of the plane of projection. If not set, an arbitrary grid-aligned 
  ``north_vector`` is chosen. Ignored in the case where a particular axis (e.g., "x", "y", or "z") is 
  explicitly specified.
* ``smooth_positions`` (optional): Apply a gaussian smoothing operation to the sky positions 
  of the events. This may be useful when the binned events appear blocky due to their uniform
  distribution within simulation cells. However, this will move the events away from their 
  originating position on the sky, and so may distort surface brightness profiles and/or 
  spectra. Should probably only be used for visualization purposes. Supply a float here to 
  smooth with a standard deviation with this fraction of the cell or particle size.
  Default: None
* ``prng`` (optional): An integer seed, pseudo-random number generator, :class:`~numpy.random.RandomState` 
  object, or :mod:`~numpy.random` (the default). Use this if you have a reason to generate the same 
  set of random numbers, such as for a test. 

Assuming one then has a :class:`~pyxsim.photon_list.PhotonList` ``photons``, example invocations could look
like this:

A simple projection along an axis:

.. code-block:: python

    events = photons.project_photons("z", (30.0, 45.0))
        
An off-axis projection:

.. code-block:: python

    events = photons.project_photons([0.1, -0.3, 0.5], (30.0, 45.0))

An on-axis projection with absorption:

.. code-block:: python

    events = photons.project_photons("y", (12.0, -30.0), absorb_model="tbabs", nH=0.01)

An off-axis projection with a ``north_vector``, without Doppler velocity shifting, 
and a specific random number generator:

.. code-block:: python
    
    events = photons.project_photons([0.1, -0.3, 0.5], (12.0, -30.0), no_shifting=True, 
                                     north_vector=[1.0,0.0,0.0], prng=34)

.. note::

    Unlike the ``photon_simulator`` analysis module in yt, the ability to convolve 
    the event energies using an ARF and RMF has been taken out of this step entirely 
    and moved into a new instrument simulator step. See :ref:`instruments` for details. 
    
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
method. Since the :class:`~pyxsim.event_list.EventList` does not have an
intrinsic binning, we need to provide a field of view ``fov`` and a resolution
``nx``:

.. code-block:: python

    fov = (10.0, "arcmin") # the field of view / width of the image
    nx = 256 # The resolution of the image on a side
    events.write_fits_file("cluster_events.fits", fov, nx, overwrite=True)
    
The ``overwrite`` keyword argument is used to allow (or prevent) overwrites of 
files if they already exist. To read previously stored events back from disk, 
use the :meth:`~pyxsim.event_list.EventList.from_fits_file` method:

.. code-block:: python

    events = EventList.from_fits_file("cluster_events.fits")

.. _simput:

SIMPUT
++++++

An :class:`~pyxsim.event_list.EventList` can be exported to the SIMPUT file format for
reading in by other packages that simulate particular instruments, such as
`SOXS <http://hea-www.cfa.harvard.edu/~jzuhone/soxs>`_, 
`MARX <http://space.mit.edu/ASC/MARX/>`_, or `SIMX <http://hea-www.cfa.harvard.edu/simx/>`_
(see also :ref:`instruments`). This is done by calling the 
:meth:`~pyxsim.event_list.EventList.write_simput_file` method:

.. code-block:: python

    events.write_simput_file("my_great_events", overwrite=False, emin=0.1, emax=9.0)

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

Saving Derived Products from Event Lists
----------------------------------------

:class:`~pyxsim.event_list.EventList` instances can produce binned images and spectra
from their events. Both products are written in FITS format.

Images
++++++

To produce a binned image, call the :meth:`~pyxsim.event_list.EventList.write_fits_image`
method:

.. code-block:: python

    fov = (20.0, "arcmin") # the field of view / width of the image
    nx = 1024 # The resolution of the image on a side
    events.write_fits_image("myimage.fits", fov, nx, overwrite=True, 
                            emin=0.5, emax=7.0)

which writes an image binned using the ``fov`` (width in angle) and ``nx`` (resolution)
parameters to the file ``"myimage.fits"``. Set ``overwrite=True`` if the file is already 
there and you want to overwrite it. The ``emin`` and ``emax`` parameters control the 
energy range of the events which will be included in the image (default is to include 
all of the events).

Spectra
+++++++

To produce a spectrum binned on energy, call :meth:`~pyxsim.event_list.EventList.write_spectrum`. 

.. code-block:: python

    specfile = "myspec.fits" # filename to write to
    emin = 0.1 # minimum energy of spectrum
    emax = 10.0 # maximum energy of spectrum
    nchan = 2000 # number of bins in spectrum
    events.write_spectrum(specfile, emin, emax, nchan, overwrite=False)

This bins the unconvolved event energies using the ``emin``, ``emax``, and ``nchan`` 
arguments into a histogram which will be written to the file as a spectrum. As usual, 
the ``overwrite`` argument determines whether or not a file can be overwritten. 
