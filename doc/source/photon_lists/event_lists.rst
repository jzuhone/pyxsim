.. _event-lists:

Event Lists
===========

Event lists are the main data products of pyXSIM, since they represent the
synthetic observations that may be compared to and analyzed in the same fashion
as real observed data. Event lists can be produced from photon list datasets by
projecting along a particular axis. They can also be used to produce other
useful data products.

Creating a New Event List by Projecting from a Photon List
----------------------------------------------------------

An event list is first created from a photon list dataset using the
:func:`~pyxsim.photon_list.project_photons` function. This method takes the
photon list dataset, a line of sight direction, and a central position on the
sky, as well as a number of additional keyword arguments to customize the
resulting projection. The events are written to disk in the same way as
photon lists are. The arguments are:

* ``photon_prefix``: The prefix of the filename(s) containing the photon list.
  If run in serial, the filename will be "{photon_prefix}.h5", if run in
  parallel, the filenames will be "{photon_prefix}.{mpi_rank}.h5".
* ``event_prefix``: The prefix of the filename(s) which will be written to
  contain the event list. If run in serial, the filename will be
  ``"{event_prefix}.h5"``, if run in parallel, the filename will be
  ``"{event_prefix}.{mpi_rank}.h5"``.
* ``normal``: The line of sight direction to project along. Accepts either a
  coordinate axis (``"x"``, ``"y"``, or ``"z"``), or a three-vector for an
  off-axis projection, e.g. ``[1.0, -0.3, 0.24]``.
* ``sky_center``: Central RA, Dec of the events in degrees.
* ``absorb_model`` (optional): A string representing a model for foreground
  galactic absorption. The two models included in pyXSIM for absorption are:
  ``"wabs"``, `Wisconsin (Morrison and McCammon; ApJ 270, 119) <http://adsabs.harvard.edu/abs/1983ApJ...270..119M>`_,
  and ``"tbabs"``, `Tuebingen-Boulder (Wilms, J., Allen, A., & McCray, R. 2000, ApJ, 542, 914) <http://adsabs.harvard.edu/abs/2000ApJ...542..914W>`_.
  The default is no absorption--if an absorption model is chosen, the ``nH``
  parameter must also be set.
* ``nH`` (optional): The foreground galactic column density in units of
  10\ :sup:`22` atoms cm :sup:`-2`, for use when one is applying
  foreground galactic absorption.
* ``abund_table`` (optional): The abundance table to be used for abundances in the
  TBabs absorption model. Default is set in the SOXS configuration file, the default
  for which is ``"angr"``. Other options are ``"angr"``, ``"aspl"``, ``"lodd"``,
  ``"feld"``, ``"wilm"``, and ``"cl17.03"``. For the definitions of these, see
  :ref:`solar-abund-tables`.
* ``no_shifting`` (optional): If set to True, the photon energies will not be
  velocity Doppler shifted. Default False.
* ``north_vector`` (optional): A vector defining the "up" direction, e.g.
  ``[0.0, 1.0, 0.0]``. This option sets the orientation of the plane of
  projection. If not set, an arbitrary grid-aligned
  ``north_vector`` is chosen. Ignored in the case where a particular axis (e.g.,
  "x", "y", or "z") is explicitly specified.
* ``sigma_pos`` (optional): Apply a gaussian smoothing operation to the sky
  positions of the events. This may be useful when the binned events appear
  blocky due to their uniform distribution within simulation cells. However,
  this will move the events away from their originating position on the sky,
  and so may distort surface brightness profiles and/or spectra. Should probably
  only be used for visualization purposes. Supply a float here to smooth with a
  standard deviation with this fraction of the cell or particle size.
  Default: None
* ``flat_sky`` (optional): If ``True``, we assume that the sky is "flat" and
  RA, Dec positions are computed using simple linear offsets, Default: ``False``.
* ``save_los`` (optional): If ``True``, save the line-of-sight positions along
  the projection axis in units of kpc to the events list. Default: ``False``.
* ``prng`` (optional): An integer seed, pseudo-random number generator,
  :class:`~numpy.random.RandomState` object, or :mod:`~numpy.random` (the
  default). Use this if you have a reason to generate the same set of random
  numbers, such as for a test.

This function returns the number of events generated. Example invocations could
look like this:

A simple projection along an axis:

.. code-block:: python

    n_events = pyxsim.project_photons("my_photons", "my_events", "z",
                                      (30.0, 45.0))

An off-axis projection:

.. code-block:: python

    n_events = pyxsim.project_photons("my_photons", "my_events",
                                      [0.1, -0.3, 0.5], (30.0, 45.0))

An on-axis projection with absorption:

.. code-block:: python

    n_events = pyxsim.project_photons("my_photons", "my_events", "y",
                                      (12.0, -30.0), absorb_model="tbabs",
                                      nH=0.01)

An off-axis projection with a ``north_vector``, without Doppler velocity
shifting, and a specific random number generator:

.. code-block:: python

    n_events = pyxsim.project_photons("my_photons", "my_events",
                                      [0.1, -0.3, 0.5], (12.0, -30.0),
                                      no_shifting=True,
                                      north_vector=[1.0,0.0,0.0], prng=34)


Reading Event Lists from Disk
-----------------------------

Event lists are written to disk by :func:`~pyxsim.photon_list.project_photons`,
and can be read back in using the :class:`~pyxsim.event_list.EventList` class.
This class facilitates various tasks for converting events to other formats.

To read in an event list, simply provide the filename if it is a single file:

.. code-block:: python

    events = pyxsim.EventList("my_events.h5")

Or, if the filenames are split into multiple numbered files, choose the
first one:

.. code-block:: python

    events = pyxsim.EventList("my_events.0000.h5")

the others will be found automatically, as the total list of files is stored in
the first one.

For event list files created previous to pyXSIM version 4.3.0, an event list
split up into multiple files should be loaded up in one of two ways. Either
the full list of files can be provided:

.. code-block:: python

    events = pyxsim.EventList(["my_events.0001.h5","my_events.0002.h5","my_events.0003.h5"])

or a regular expression which can be used to infer the filenames:

.. code-block:: python

    events = pyxsim.EventList("my_events*.h5")

The parameters used in the run to produce the event list are stored in a
``parameters`` dictionary:

.. code-block:: python

    print(events.parameters)

.. code-block:: pycon

    {'absoption_model': 'wabs',
     'abund_table': 'angr',
     'area': 500.0,
     'exp_time': 100000.0,
     'flat_sky': 0,
     'kernel': 'top_hat',
     'nH': 0.02,
     'no_shifting': 0,
     'normal': 'x',
     'observer': 'external',
     'sky_center': array([30., 45.])}

and other pertinent information used in the production of the event list can
be found in the attached ``info`` dictionary:

.. code-block:: python

    print(events.info)

.. code-block:: pycon

    {'photon_file': 'plaw_photons.h5',
     'pyxsim_version': '4.1b1.dev29+g1c09873.d20221228',
     'soxs_version': 'soxs-4.2.2.dev22+gd56e1b4',
     'yt_version': '4.2.dev0'}

If this event list file has originated from merged event lists, then there
will be multiple instances of each piece of information, numbered by the
file, e.g. ``"soxs_version_0"``, ``"soxs_version_1"``, and so on. The original
files used to make the merge will be stored in the key ``"original_files"``.

.. _simput:

SIMPUT
++++++

An :class:`~pyxsim.event_list.EventList` can be exported to the SIMPUT file
format for reading in by other packages that simulate particular instruments,
such as `SOXS <http://hea-www.cfa.harvard.edu/soxs>`_,
`MARX <http://space.mit.edu/ASC/MARX/>`_, or
`SIMX <http://hea-www.cfa.harvard.edu/simx/>`_
(see also :ref:`instruments`). This is done by calling the
:meth:`~pyxsim.event_list.EventList.write_simput_file` method:

.. code-block:: python

    events.write_simput_file("my_great_events", overwrite=False,
                             emin=0.1, emax=9.0)

where the first argument is the prefix for the files that will be created (the
SIMPUT file and a photon list sidecar file), and the other optional argument
controls whether or not an existing file will be overwritten. Currently, SIMPUT
files are used for export only; they cannot be used to read events back into
pyXSIM.

Images
++++++

To produce a binned image, call the
:meth:`~pyxsim.event_list.EventList.write_fits_image` method:

.. code-block:: python

    fov = (20.0, "arcmin") # the field of view / width of the image
    nx = 1024 # The resolution of the image on a side
    events.write_fits_image("myimage.fits", fov, nx, overwrite=True,
                            emin=0.5, emax=7.0)

which writes an image binned using the ``fov`` (width in angle) and ``nx``
(resolution) parameters to the file ``"myimage.fits"``. Set ``overwrite=True``
if the file is already there and you want to overwrite it. The ``emin`` and
``emax`` parameters control the energy range of the events which will be
included in the image (default is to include all of the events).

Spectra
+++++++

To produce a spectrum binned on energy, call
:meth:`~pyxsim.event_list.EventList.write_spectrum`.

.. code-block:: python

    specfile = "myspec.fits" # filename to write to
    emin = 0.1 # minimum energy of spectrum
    emax = 10.0 # maximum energy of spectrum
    nchan = 2000 # number of bins in spectrum
    events.write_spectrum(specfile, emin, emax, nchan, overwrite=False)

This bins the unconvolved event energies using the ``emin``, ``emax``, and
``nchan`` arguments into a histogram which will be written to the file as a
spectrum. As usual, the ``overwrite`` argument determines whether or not a file
can be overwritten.

Merging Event Lists
-------------------

Event lists which have been written to files can be merged together, using the
:func:`~pyxsim.utils.merge_files` function. This may be useful if you generate events from
different sources or source types that are co-located on the sky.

:func:`~pyxsim.utils.merge_files` takes a list of input filenames, and an output filename.
The optional keyword arguments are ``overwrite``, which decides whether or not an existing file
will be overwritten, and ``add_exposure_times`` decides whether or not the final file will
have an exposure time of the sum of the times in the separate files or that of the longest
exposure time between the files.

.. code-block:: python

    from pyxsim import merge_files
    merge_files(["events_0.h5","events_1.h5","events_3.h5"], "events.h5",
                overwrite=True, add_exposure_times=True)
