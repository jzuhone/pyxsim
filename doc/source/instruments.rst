.. _instruments:

Convolving Events with Instrumental Responses
=============================================

Built-In Convolutions
---------------------

pyXSIM provides the ability to perform "quick-and-dirty" convolutions in energy
and spatial position. The accuracy of these methods is limited, since they

Convolving with an ARF
++++++++++++++++++++++

An energy-dependent effective area can be taken into account in the creation of
the :class:`~pyxsim.event_list.EventList` using :meth:`~pyxsim.photon_list.PhotonList.project_photons`. 
Simply create a :class:`~pyxsim.responses.AuxiliaryResponseFile` object, and 
pass it in as the ``area_new`` parameter to :meth:`~pyxsim.photon_list.PhotonList.project_photons`:

.. code-block:: python

.. note::

    If you want to export simulated events to an instrument simulator, DO NOT
    perform any of these convolution methods within pyXSIM. This will
    be the job of the instrument simulator itself to do.  

MARX
----

As of MARX version 5.3, simulated events can be read in via ``SIMPUT`` files. 


SIMX
----

Events can be inputted into ``SIMX`` by use of a ``SIMPUT`` file. 