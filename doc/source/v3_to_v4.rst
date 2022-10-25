.. _v3_to_v4:

Important API Changes Between pyXSIM Version 3 and Version 4
============================================================

This section documents the most important API changes between
pyXSIM version 3 and 4, with code examples showing how to update 
your code for version 4. 

API Changes to Thermal Source Models
------------------------------------

The class ``ThermalSourceModel`` in version 3 of pyXSIM has now been
subclassed into three separate classes, 
:class:`~pyxsim.source_models.thermal_sources.CIESourceModel`,
:class:`~pyxsim.source_models.thermal_sources.NEISourceModel`, and
:class:`~pyxsim.source_models.thermal_sources.IGMSourceModel`.
``ThermalSourceModel`` itself should no longer be used as a standalone
class. In order to replicate the functionality of the previously existing
``ThermalSourceModel`` class for CIE spectra, change (for example):

.. code-block:: python

    import pyxsim
    
    source_model = pyxsim.ThermalSourceModel("apec", 0.1, 10.0, 5000, 0.3)

to:

.. code-block:: python

    import pyxsim
    
    source_model = pyxsim.CIESourceModel("apec", 0.1, 10.0, 5000, 0.3)

and for NEI sources, change:

.. code-block:: python

    import pyxsim
    
    source_model = pyxsim.ThermalSourceModel("apec", 0.1, 10.0, 5000, 
                                             ("gas", "metallicity"), var_elem=var_elem,
                                             nei=True)
                                             
to:

.. code-block:: python

    import pyxsim
    
    source_model = pyxsim.NEISourceModel(0.1, 10.0, 5000, var_elem)

See :ref:`thermal-sources` for more details.