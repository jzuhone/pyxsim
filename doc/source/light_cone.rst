.. _light-cone:

Generating Light Cone Simulations of X-rays
===========================================

.. warning::

    This feature is currently in beta, and only works for Enzo cosmological
    simulations.

First, one needs to create the :class:`~pyxsim.light_cone.XrayLightCone` object:

.. code-block:: python

    lc = pyxsim.XrayLightCone('Enzo_64/64.param', 'Enzo', 0.0, 0.9)
    
This simply takes a simulation

After the :class:`~pyxsim.light_cone.XrayLightCone` has been created, we have to 
implement a source model to 

.. code-block:: python

    source_model = pyxsim.ThermalSourceModel("apec", 0.1, 10.0, 1000)

.. code-block:: python

    exp_time = 50000.0 # exposure time in seconds
    area = 25000.0 # collecting area in cm**2
    fov = (2.0, "deg") # field of view
    sky_center = (30.0, 45.0) # sky center in degrees
    events = lc.generate_events(area, exp_time, fov, source_model, 
                                sky_center, absorb_model="wabs", 
                                nH=0.02, smooth_positions=0.5)