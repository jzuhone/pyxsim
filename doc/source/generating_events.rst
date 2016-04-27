.. _event-lists:

Generating Event Lists
======================

Projecting from a Photon List
-----------------------------

Saving/Reading Events to/from Disk
----------------------------------

Any :class:`~pyxsim.event_list.EventList` instance may be saved to disk 
by calling the :meth:`~pyxsim.event_list.EventList.write_h5_file` method:

.. code-block:: python
    
    events.write_h5_file("cluster_events.h5")
    
This writes the event data and associated parameters to disk in the convenient HDF5 
file format. To read previously stored events back from disk, use the 
:meth:`~pyxsim.event_list.EventList.from_h5_file` method:

.. code-block:: python

    events = EventList.from_h5_file("cluster_events.h5")

Adding Background and Point Source Events
-----------------------------------------

Manipulating Event Lists
------------------------
