.. _instruments:

Convolving Events with Instrumental Responses
=============================================

Built-In Instrument Simulators
------------------------------

pyXSIM provides the ability to perform "quick-and-dirty" convolutions with approximate
representations of real X-ray instruments. The accuracy of these representations is 
limited, since they do not take into account other effects, such as vignetting, the 
actual layout of the CCDs, and position-dependent effective area and spectral response.
However, they should provide approximate representations of X-ray observations that will
be sufficient for most purposes. 

.. note::

    If you want to export simulated events to an external instrument simulator
    such as MARX or SIMX, DO NOT use any of these built-in instrument simulators
    within pyXSIM. This will be the job of the instrument simulator itself to do.

Using a Built-In Instrument Simulator
+++++++++++++++++++++++++++++++++++++

pyXSIM comes with the following built-in instrument simulators:

* ``ACIS_I``: ACIS-I Cycle 17 on-axis ARF and RMF, with 0.492" pixels 
* ``ACIS_S``: ACIS-S Cycle 17 on-axis ARF and RMF, with 0.492" pixels 
* ``AstroH_SXS``:
* ``XRS_Imager``:
* ``XRS_Calorimeter``:

Creating Your Own Built-In Instrument Simulator
+++++++++++++++++++++++++++++++++++++++++++++++

Re-Binning the Events
+++++++++++++++++++++

Convolving with a Spatial PSF
+++++++++++++++++++++++++++++

Convolving with an ARF
++++++++++++++++++++++

An energy-dependent effective area can be taken into account in the creation of
the :class:`~pyxsim.event_list.EventList` using :meth:`~pyxsim.photon_list.PhotonList.project_photons`. 
Simply create a :class:`~pyxsim.responses.AuxiliaryResponseFile` object, and 
pass it in as the ``area_new`` parameter to :meth:`~pyxsim.photon_list.PhotonList.project_photons`:

.. code-block:: python

    import pyxsim
    arf = pyxsim.AuxiliaryResponseFile("aciss_aimpt_cy17.arf")
    ...
    events = photons.project_photons("z", area_new=arf, exp_time_new=(50.,"ks"), 
                                     absorb_model=tbabs_model, sky_center=(45.,30.))

some ARFs, particularly those from Japanese missions, are not "normalized". In order to
use these, the corresponding RMF should be provided as well:

.. code-block:: python

    import pyxsim
    arf = pyxsim.AuxiliaryResponseFile("sxt-s_120210_ts02um_intallpxl.arf", 
                                       rmffile="ah_sxs_5ev_basefilt_20100712.rmf")
    ...
    events = photons.project_photons("z", area_new=arf, exp_time_new=(50.,"ks"), 
                                     absorb_model=tbabs_model, sky_center=(45.,30.))

.. note::
    
    Regardless of the origin of the ARF, the effective area curve therein is assumed
    to be spatially constant. 
    
Convolving with an RMF
++++++++++++++++++++++
    
Producing More Realistic Observations Using External Packages
-------------------------------------------------------------


MARX
++++

As of MARX version 5.3, simulated events can be read in via ``SIMPUT`` files. 


SIMX
++++

Events can be inputted into ``SIMX`` by use of a ``SIMPUT`` file. 