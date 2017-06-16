.. _absorb-models:

Absorption Models
=================

Foreground galactic absorption can be applied during the creation of events. The two 
models included in pyXSIM for absorption are:

* ``"wabs"``: Wisconsin (Morrison and McCammon; ApJ 270, 119) absorption model
* ``"tbabs"``: Tuebingen-Boulder (Wilms, J., Allen, A., & McCray, R. 2000, ApJ, 542, 914)
  absorption model

An absorption model may be specified in the call to 
:meth:`~pyxsim.photon_list.PhotonList.project_photons`. The other required parameter is
``nH``, the neutral hydrogen column density in units of :math:`10^{22} \rm{atoms} \rm{cm}^{-2}`:

.. code-block:: python

    events = photons.project_photons("y", (12.0, -30.0), absorb_model="wabs", nH=0.04)

If one has their own absorption model that they would like to use, you can create a new
absorption model class based on pyXSIM's :class:`~pyxsim.spectral_models.TableAbsorbModel` 
class and supply it to :meth:`~pyxsim.photon_list.PhotonList.project_photons` instead. 
In this case, a table of the cross sections is required, to be stored in an HDF5 file. 
This file must have two top-level datasets:

* ``"energy"``: A 1-D array of M+1 energies in units of keV, where M is the number of bins
* ``"cross_section"``: A 1-D array of M cross-sections in units of :math:`\rm{cm}^2`, where M
  is the number of bins

An example file is provided with pyXSIM, ``tbabs_table.h5``, containing a table of the 
`tbabs <http://pulsar.sternwarte.uni-erlangen.de/wilms/research/tbabs/>`_ cross sections,
generated from XSPEC. 

A new absorption model class can be created in this way:

.. code-block:: python

    from pyxsim import TableAbsorbModel

    class MyAbsorbModel(TableAbsorbModel):
        def __init__(self, nH):
            super(TBabsModel, self).__init__("absorb_table.h5", nH)
    
Then, the name of the new class can be supplied to 
:meth:`~pyxsim.photon_list.PhotonList.project_photons`:

.. code-block:: python

    events = photons.project_photons("y", (12.0, -30.0), 
                                     absorb_model=MyAbsorbModel, nH=0.04)

