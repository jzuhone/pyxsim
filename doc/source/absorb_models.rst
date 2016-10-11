.. _absorb-models:

Absorption Models
=================

Foreground galactic absorption can be applied during the creation of events. An absorption 
model must be created using one of the classes below and passed to 
:meth:`~pyxsim.photon_list.PhotonList.project_photons`, 
:meth:`~pyxsim.event_list.EventList.add_background`,
or :meth:`~pyxsim.event_list.EventList.add_point_sources`. The different ways of constructing
absorption models are:

:class:`~pyxsim.spectral_models.WabsModel` generates a Wisconsin (Morrison and McCammon; 
ApJ 270, 119) absorption model ("wabs"):

.. code-block:: python

    N_H = 0.1 # galactic column density in units of 10^{22} cm^{-2}
    abs_model = pyxsim.WabsModel(N_H)

:class:`~pyxsim.spectral_models.TBabsModel` generates a Tuebingen-Boulder (Wilms, J., 
Allen, A., & McCray, R. 2000, ApJ, 542, 914) absorption model ("TBabs"):

.. code-block:: python

    N_H = 0.1 # galactic column density in units of 10^{22} cm^{-2}
    abs_model = pyxsim.TBabsModel(N_H)

:class:`~pyxsim.spectral_models.XSpecAbsorbModel` generates an absorption spectrum from 
an XSpec model (must have PyXspec installed):

.. code-block:: python

    model = "wabs" # or "phabs", or "TBabs", etc.
    N_H = 0.1 # galactic column density in units of 10^{22} cm^{-2}
    abs_model = pyxsim.XSpecAbsorbModel(model, N_H)

:class:`~pyxsim.spectral_models.TableAbsorbModel` generates an absorption spectrum from 
an HDF5-based table of energy and cross section:

.. code-block:: python

    filename = "tbabs_table.h5" # file containing the table
    N_H = 0.1 # galactic column density in units of 10^{22} cm^{-2}
    abs_model = pyxsim.XSpecAbsorbModel(filename, N_H)
    
The HDF5 file must have two top-level datasets:
 
* ``"energy"``: A 1-D array of M+1 energies in units of keV, where M is the number of bins
* ``"cross_section"``: A 1-D array of M cross-sections in units of :math:`\rm{cm}^2`, where M
  is the number of bins
  
An example file is provided with pyXSIM, ``tbabs_table.h5``, containing a table of the 
`tbabs <http://pulsar.sternwarte.uni-erlangen.de/wilms/research/tbabs/>`_ cross sections,
generated from XSPEC. 

