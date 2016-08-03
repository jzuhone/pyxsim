.. _spectral-models:

Special-Purpose Spectral Models
===============================

pyXSIM provides special classes to handle specific types of spectral models that are not reducible
to simply analytical models. These classes read the models either directly from tables or process
them using PyXspec. Currently, this includes models representing emission from thermal plasmas and 
galactic foreground absorption models.

Thermal Models
--------------

For thermal plasmas, pyXSIM provides classes for modeling emission that is dependent on temperature
and metallicity, e.g. :math:`\Lambda(T,Z)`. These can serve as input to the 
:class:`~pyxsim.source_models.ThermalSourceModel` class. 

:class:`~pyxsim.spectral_models.XSpecThermalModel` generates a thermal emission spectrum
from a model known to XSPEC, using PyXspec as a backend:

.. code-block:: python

    model = "apec" # or "mekal", "bapec", etc.
    emin = 0.01 # The minimum energy of the spectrum in keV
    emax = 20.0 # The maximum energy of the spectrum in keV
    nchan = 10000 # The number of spectral channels
    spec_model = pyxsim.XSpecThermalModel(model, emin, emax, nchan, thermal_broad=True)

They keyword argument ``thermal_broad`` should be set to ``True`` or ``False`` depending on
whether or not you want the spectral lines thermally broadened. 

:class:`~pyxsim.spectral_models.TableApecModel` generates a thermal emission spectrum
from the APEC plasma emission tables available from `AtomDB <http://www.atomdb.org>`_:

.. code-block:: python

    emin = 0.01 # The minimum energy of the spectrum in keV
    emax = 20.0 # The maximum energy of the spectrum in keV
    nchan = 10000 # The number of spectral channels
    apec_vers = "3.0.3" # The version identifier string for the APEC files. Default: "2.0.2"
    apec_root = "/Users/jzuhone/atomdb" # The directory where the APEC model files are stored.
                                        # Optional, the native pyxsim files will be used if
                                        # a location is not provided.
    spec_model = pyxsim.TableApecModel(emin, emax, nchan, apec_root=apec_root,
                                       apec_vers=apec_vers, thermal_broad=False)

You will need to set up one of these two models in your script and pass it as the first argument to
:class:`~pyxsim.source_models.ThermalSourceModel` (see :ref:`thermal-sources` for details).

Absorption Models
-----------------

Foreground galactic absorption is applied during the creation of events, whether within 
:meth:`~pyxsim.photon_list.PhotonList.project_photons`, :meth:`~pyxsim.event_list.EventList.add_background`,
or :meth:`~pyxsim.event_list.EventList.add_point_sources`. There are two classes in pyXSIM
for generating foreground absorption spectra. 

:class:`~pyxsim.spectral_models.XSpecAbsorbModel` generates an absorption spectrum from 
an XSpec model (must have PyXspec installed):

.. code-block:: python

    model = "wabs" # or "phabs", or "TBabs", etc.
    N_H = 0.1 # galactic column density in 10^{22} cm^{-2}
    abs_model = pyxsim.XSpecAbsorbModel(model, N_H)

:class:`~pyxsim.spectral_models.TableAbsorbModel` generates an absorption spectrum from 
an HDF5-based table of energy and cross section:

.. code-block:: python

    filename = "tbabs_table.h5" # file containing the table
    N_H = 0.1 # galactic column density in 10^{22} cm^{-2}
    abs_model = pyxsim.XSpecAbsorbModel(model, N_H)
    
The HDF5 file must have two top-level datasets:
 
* ``"energy"``: A 1-D array of M+1 energies in units of keV, where M is the number of bins
* ``"cross_section"``: A 1-D array of M cross-sections in units of :math:`\rm{cm}^2`, where M
  is the number of bins
  
An example file is provided with pyXSIM, ``tbabs_table.h5``, containing a table of the 
`tbabs <http://pulsar.sternwarte.uni-erlangen.de/wilms/research/tbabs/>`_ cross sections,
generated from XSPEC. 

