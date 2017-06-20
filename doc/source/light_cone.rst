.. _light-cone:

Generating Light Cone Simulations of X-rays
===========================================

.. warning::

    This feature is currently in beta, and only works for Enzo cosmological
    simulations. In the future, it will be expanded to other types of 
    datasets. If you'd like to see support for your type of dataset, 
    `get in contact <faq.html>`_!

First, one needs to create the :class:`~pyxsim.light_cone.XrayLightCone` object:

.. code-block:: python

    lc = pyxsim.XrayLightCone('Enzo_64/64.param', 'Enzo', 0.0, 0.9)
    
This simply takes a simulation

After the :class:`~pyxsim.light_cone.XrayLightCone` has been created, we have to 
implement a source model to determine how the photons will be generated from the
source properties (as usual). In this case, we'll simply use the 
:class:`~pyxsim.source_models.ThermalSourceModel`:

.. code-block:: python

    source_model = pyxsim.ThermalSourceModel("apec", 0.1, 10.0, 1000)

Now, we are ready to generate our simulated events using 
:meth:`~pyxsim.light_cone.XrayLightCone.generate_events`. Since by definition 
a light cone is a projection, in this case there is no intermediate step of 
creating a :class:`~pyxsim.photon_list.PhotonList` first--the output is an
:class:`~pyxsim.event_list.EventList`. We must specify an exposure time, 
collecting area, field of view in units of angle, the source model, and the 
center of the field of view in (RA, Dec). The 
:meth:`~pyxsim.light_cone.XrayLightCone.generate_events` 
method also takes a number of the same optional parameters as 
:class:`~pyxsim.photon_list.PhotonList.project_photons`, so we'll also absorb
the events with the ``wabs`` model and smooth the positions just a bit for
visualization purposes.

.. code-block:: python

    exp_time = 50000.0 # exposure time in seconds
    area = 25000.0 # collecting area in cm**2
    fov = (2.0, "deg") # field of view
    sky_center = (30.0, 45.0) # sky center in degrees
    events = lc.generate_events(area, exp_time, fov, source_model, 
                                sky_center, absorb_model="wabs", 
                                nH=0.02, smooth_positions=0.5)

If we make an image of this :class:`~pyxsim.event_list.EventList`, it looks
like this: