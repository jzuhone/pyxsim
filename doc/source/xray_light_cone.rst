.. _xray-light-cone:

Generating Light Cone Simulations of X-rays
===========================================

.. warning::

    This feature is currently in beta, and only works for Enzo cosmological
    simulations.

First, one needs to create the :class:`~pyxsim.light_cone.XrayLightCone` object:

.. code-block:: python

    lc = pyxsim.XrayLightCone('Enzo_64/64.param', 'Enzo', 0.0, 0.9)
    
This simply takes a simulation

After the :class:`~pyxsim.light_cone.XrayLightCone` has been created, 
