.. _absorb-models:

Absorption Models
=================

Foreground galactic absorption can be applied during the creation of events. 
The two models included in pyXSIM for absorption are:

* ``"wabs"``: `Wisconsin (Morrison and McCammon; ApJ 270, 119) <http://adsabs.harvard.edu/abs/1983ApJ...270..119M>` 
  absorption model
* ``"tbabs"``: `Tuebingen-Boulder (Wilms, J., Allen, A., & McCray, R. 2000, ApJ, 542, 914) <http://adsabs.harvard.edu/abs/2000ApJ...542..914W>`_
  absorption model

An absorption model may be specified in the call to 
:meth:`~pyxsim.photon_list.project_photons`. The other required 
parameter is ``nH``, the neutral hydrogen column density in units of 
:math:`10^{22} \rm{atoms}~\rm{cm}^{-2}`:

.. code-block:: python

    n_events = pyxsim.project_photons("my_photons", "my_events", "y", 
                                      (12.0, -30.0), absorb_model="wabs", 
                                      nH=0.04)

If one has their own absorption model that they would like to use, you can 
create a new absorption model class based on pyXSIM's 
:class:`~pyxsim.spectral_models.TableAbsorbModel` class and supply it to 
:meth:`~pyxsim.photon_list.project_photons` instead. In this case, a table of 
the cross sections is required, in the form of two NumPy arrays, which should 
be:

* ``energy``: A 1-D array of M energies in units of keV, where M is the number 
  of bins
* ``cross_section``: A 1-D array of M cross-sections in units of 
  :math:`\rm{cm}^2`, where M is the number of bins

A new absorption model class can be created in this way:

.. code-block:: python

    from pyxsim import AbsorptionModel

    # define energy and cross_section as NumPy arrays somewhere here
    energy = ...
    cross_section = ...
 
    class MyAbsorbModel(AbsorptionModel):
        def __init__(self, nH):
            super(AbsorptionModel, self).__init__(nH, energy, cross_section)
    
Then, the name of the new class can be supplied to 
:meth:`~pyxsim.photon_list.project_photons`:

.. code-block:: python

    n_events = pyxsim.project_photons("my_photons", "my_events", "y", 
                                      (12.0, -30.0), 
                                      absorb_model=MyAbsorbModel, nH=0.04)

