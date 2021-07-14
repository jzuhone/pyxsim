.. _v1_to_v2:

Important API Changes Between pyXSIM Version 2 and Version 3
============================================================

This section documents the most important API changes between
pyXSIM version 2 and 3, with code examples showing how to update 
your code for version 3. 

Metallicity in Thermal Sources
------------------------------

In the :class:`~pyxsim.source_models.ThermalSourceModel` class, the ``Zmet``
keyword argument has been changed to a required argument. Whereas previously 
one would not have to set ``Zmet`` at all, and it would default to 
0.3 :math:`Z_\odot`, it now must be set. 

Old way:

.. code-block:: python

    # spatially constant metallicity
    source_model = pyxsim.ThermalSourceModel("apec", 0.1, 11.0, 1000.0, Zmet=0.4)
    
    # metallicity field
    source_model = pyxsim.ThermalSourceModel("apec", 0.1, 11.0, 1000.0, 
                                             Zmet=("gas","metallicity"))

New way:

.. code-block:: python

    # spatially constant metallicity
    source_model = pyxsim.ThermalSourceModel("apec", 0.1, 11.0, 1000.0, 0.4)
    
    # metallicity field
    source_model = pyxsim.ThermalSourceModel("apec", 0.1, 11.0, 1000.0, 
                                             ("gas","metallicity"))

Creating Photon Lists
---------------------

Photon lists are no longer created as a ``PhotonList`` class, but are written
directly to disk as the photons are generated. The change was made since photon 
lists can sometimes become quite large and difficult to fit into memory, even if
it does not take very long to generate them. The effect is essentially the same
as if in pyXSIM 2.x one were to generate a ``PhotonList`` and write it immediately
to disk.

Old way of creating a photon list:

.. code-block:: python
 
    photons = pyxsim.PhotonList.from_data_source(sp, redshift, area, exp_time,
                                                 source_model, center=center, 
                                                 cosmology=cosmo)
 
    photons.write_h5_file("my_photons.h5")
    
New way of creating a photon list, with the identical result:

.. code-block:: python

    n_photons, n_cells = pyxsim.make_photons("my_photons", sp, redshift, area,
                                             exp_time, source_model, 
                                             center=center, cosmology=cosmo)

All other optional arguments are the same. See :ref:`photon-lists` for more
information. 

Creating Event Lists
--------------------

Similarly, event lists (projected, redshifted, absorbed photons) are written
directly to disk as they are created, and are created from photons which are
also read from disk in chunks. This results in a similar API change:

Old way of creating an event list:

.. code-block:: python

    photons = pyxsim.PhotonList.from_file("my_photons.h5")
    events = photons.project_photons("z", (30.0, 45.0))
    events.write_h5_file("my_events.h5")
    
New way of creating an event list:

.. code-block:: python
    
    n_events = pyxsim.project_photons("my_photons", "my_events", "z", 
                                      (30.0, 45.0))
    
All optional arguments are the same. See :ref:`event-lists` for more 
information. 

Using Event Lists
-----------------

There is still a :class:`~pyxsim.event_list.EventList` class in pyXSIM 3.x. 
To create an :class:`~pyxsim.event_list.EventList` instance, the only way
to do it now is to read it from disk:

.. code-block:: python

    events = pyxsim.EventList("my_events.h5")

The way to write a SIMPUT catalog from an :class:`~pyxsim.event_list.EventList` 
has changed slightly. The old way was:

.. code-block:: python

    events.write_simput_file("my_great_events", overwrite=False)

The new way is:

.. code-block:: python

    events.write_to_simput("my_great_events", overwrite=False)

See :ref:`event-lists` for more information. 

Generating Background and Point Source Events
---------------------------------------------

The source generator functions :func:`~pyxsim.source_generators.background.make_background`
and :func:`~pyxsim.source_generators.point_sources.make_point_sources`
no longer exist in pyXSIM. To make background and point source events, 
please consult the `SOXS <https://hea-www.cfa.harvard.edu/soxs>`_ package.
