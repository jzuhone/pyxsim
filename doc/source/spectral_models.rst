.. _spectral-models:

Special-Purpose Spectral Models
===============================

pyXSIM provides special classes to handle specific types of spectral models that may be used often. 

Thermal Models
--------------

Absorption Models
-----------------

Foreground galactic absorption is applied during the creation of events, whether within 
:meth:`~pyxsim.photon_list.PhotonList.project_photons`, :meth:`~pyxsim.event_list.EventList.add_background`,
or :meth:`~pyxsim.event_list.EventList.add_point_sources`. There are two classes in pyXSIM
for generating absorption models. 

:class:`~pyxsim.spectral_models.XSpecAbsorbModel` generates an absorption model from 
an XSpec model (must have PyXspec installed):

.. code-block:: python

    model = "wabs" # or "phabs", or "tbabs", etc.
    N_H = 0.1 # galactic column density in 10^{22} cm^{-2}
    abs_model = XSpecAbsorbModel(model, N_H)

:class:`~pyxsim.spectral_models.TableAbsorbModel` generates an absorption model from 
an HDF5-based table of energy and cross section:

.. code-block:: python

    filename = "tbabs_table.h5" # file containing the table
    N_H = 0.1 # galactic column density in 10^{22} cm^{-2}
    abs_model = XSpecAbsorbModel(model, N_H)
    
The HDF5 file must have two top-level datasets:
 
* ``"energy"``: A 1-D array of M+1 energies in units of keV, where M is the number of bins
* ``"cross_section"``: A 1-D array of M cross-sections in units of :math:`\rm{cm}^2`, where M
  is the number of bins
  
An example file is provided with pyXSIM, ``tbabs_table.h5``, containing a table of the 
`tbabs <http://pulsar.sternwarte.uni-erlangen.de/wilms/research/tbabs/>`_ cross sections
generated from XSPEC. 

