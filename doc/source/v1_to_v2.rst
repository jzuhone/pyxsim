.. _v1_to_v2:

Important API Changes Between pyXSIM Version 1 and Version 2
============================================================



Projecting Photons from a Photon List
-------------------------------------

Generating Background and Point Source Events
---------------------------------------------

The :class:`~pyxsim.event_list.EventList` methods ``add_background`` and
``add_point_sources`` have been replaced by the new source generator
functions :func:`~pyxsim.source_generators.background.make_background`
and :func:`~pyxsim.source_generators.point_sources.make_point_sources`.

Instrument Simulators
---------------------

The only operations now performed by :class:`~pyxsim.instruments.InstrumentSimulator` are
convolution with the effective area curve (using the ARF) and with the response matrix
(using the RMF). No spatial PSF convolutions or rebinning operations can be applied. For
more detailed instrument simulation, users are advised to write events to SIMPUT files
and use SOXS directly. 

:class:`~pyxsim.instruments.InstrumentSimulator` instances now only take
a single argument, the :class:`~pyxsim.event_list.EventList`:

.. code-block:: python

    convolved_events = pyxsim.ACIS_I(events)
    
The object which is returned is a :class:`~pyxsim.event_list.ConvolvedEventList`
which has methods specific to dealing with spectral channels. See :ref:`convolved_events` 
for more information.

