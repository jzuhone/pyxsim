.. _v1_to_v2:

Important API Changes Between pyXSIM Version 1 and Version 2
============================================================

This section documents the most important API changes between
pyXSIM version 1 and 2, with code examples showing how to update 
your code for version 2. 

Creating Thermal Source Models
------------------------------

Old calling signature example:

.. code-block:: python

New calling signature example:

.. code-block:: python

Projecting Photons from a Photon List
-------------------------------------

Two major changes to :meth:`~pyxsim.PhotonList.project_photons`
are important to note:

* The ``sky_center`` argument is now a required argument. 
* The ``absorb_model`` argument now takes a string, ``"wabs"`` or 
  ``"tbabs"``, instead of a :class:`~pyxsim.spectral_models.AbsorptionModel` 
  instance, and ``nH`` has been added as an optional argument
  to :meth:`~pyxsim.PhotonList.project_photons`. It is no longer
  necessary to create the :class:`~pyxsim.spectral_models.AbsorptionModel` 
  on your own, as :meth:`~pyxsim.PhotonList.project_photons`
  now does this internally (however, it is still possible to pass
  the specific :class:`~pyxsim.spectral_models.AbsorptionModel` class
  itself, which allows for custom absorption models to be implemented.)  

Old calling signature example:

.. code-block:: python

    nH = 0.02
    sky_center = (30.0, 45.0)

    tbabs_model = pyxsim.TBabsModel(nH)
    
    events = photons.project_photons("z", absorb_model=tbabs_model, 
                                     sky_center=sky_center)
                                     
New calling signature example:

.. code-block:: python

    nH = 0.02
    sky_center = (30.0, 45.0)
 
    events = photons.project_photons("z", sky_center, absorb_model="tbabs", 
                                     nH=nH)
    
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

