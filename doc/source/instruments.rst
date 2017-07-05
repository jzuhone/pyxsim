.. _instruments:

Convolving Events with Instrumental Responses
=============================================

Built-In Instrument Simulators
------------------------------

pyXSIM provides the ability to perform "quick-and-dirty" convolutions with 
approximate representations of real X-ray instruments using the 
`SOXS <http://hea-www.cfa.harvard.edu/~jzuhone/soxs>`_ software package as a 
backend. The accuracy of these representations is limited, since they 
assume the following simplifications:

* The field of view and angular resolution of the observation is the same
  as the unconvolved events (e.g., whatever they are from the particular
  data object you generated the events from). No binning of the event
  positions based on actual values of the field of view and/or angular
  resolution of the telescope is performed.
* The spectral response and effective area are position-independent.
* No spatial PSF is applied. 

If you only need to represent how many X-ray counts you expect for a given
observation and perhaps make spectra, using a built-in instrument simulator 
should be sufficient for your purposes. For anything more complicated, using
one of the software packages mentioned below in :ref:`realism` is recommended.

Using a Built-In Instrument Simulator
+++++++++++++++++++++++++++++++++++++

pyXSIM comes with the following built-in instrument simulators:

* ``ACIS_I``: ACIS-I Cycle 18 on-axis ARF and RMF
* ``ACIS_S``: ACIS-S Cycle 18 on-axis ARF and RMF
* ``Hitomi_SXS``: Hitomi Soft X-ray Spectrometer ARF and RMF
* ``Athena_WFI``: Athena Wide-Field Imager ARF and RMF
* ``Athena_XIFU``: Athena X-ray Integral Field Unit ARF and RMF
* ``Lynx_Imager``: X-ray Surveyor Imager ARF and RMF
* ``Lynx_Calorimeter``: X-ray Surveyor Calorimeter ARF and RMF

When an instrument simulator is called, the following operations are applied 
to the input events, in this order:

1. Using the effective area curve from the selected ARF, events are selected 
   or rejected for observation.
2. The observed event energies are convolved with the selected RMF to produce 
   the observed energy channels. 

An ``InstrumentSimulator`` takes an :class:`~pyxsim.event_list.EventList`, performs 
the above operations, and returns a :class:`~pyxsim.event_list.ConvolvedEventList`, 
containing the events which remain after being "observed" with the instrument's effective 
area and the convolved channel energies (PI or PHA depending on the RMF). For more 
information about the convolved events, see :ref:`convolved_events`. Using an 
``InstrumentSimulator`` is very simple:

.. code-block:: python

    >>> from pyxsim import ACIS_S
    >>> new_events = ACIS_S(events)

Then, event files or spectra can be written:

.. code-block:: python

    >>> new_events.write_fits_file("acis_events.fits", overwrite=True)
    >>> new_events.write_channel_spectrum("acis_spec.pi", overwrite=True)

Designing Your Own Instrument Simulator
+++++++++++++++++++++++++++++++++++++++

If you want to design an instrument simulator yourself for use with pyXSIM,  
it is fairly simple. You need only to provide the following parameters in the
call to ``InstrumentSimulator``, in this order: 

* ``name``: The short-hand name for the instrument.
* ``arf_file``: The path to the ARF file you want to use. 
* ``rmf_file``: The path to the RMF file you want to use. 

.. warning::

    It goes without saying that the ARF and RMF you choose MUST be consistent with each other; e.g., 
    have the same energy binning, etc.
    
Here we show an example of how to set up an instrument simulator, using the parameters that pyXSIM
uses for the built-in ACIS-S simulator:

.. code-block:: python

    from pyxsim import InstrumentSimulator

    ACIS_S = InstrumentSimulator("acis-s", "aciss_aimpt_cy18.arf",
                                 "aciss_aimpt_cy18.rmf")

.. _realism:

Producing More Realistic Observations Using External Packages
-------------------------------------------------------------

If you want to produce a more realistic simulation of a particular instrumental configuration,
pyXSIM provides options for exporting its event lists to external packages. For 
`MARX <http://space.mit.edu/ASC/MARX/>`_ and `SIMX <http://hea-www.cfa.harvard.edu/simx/>`_, or
for more fine-tuned use of `SOXS <http://hea-www.cfa.harvard.edu/~jzuhone/soxs>`_, one can use 
SIMPUT files. 

MARX
++++

The MARX version needs to be at least 5.3.1. To use SIMPUT with MARX, one only needs to 
change the following lines in the ``marx.par`` file:

.. code-block:: bash

    # Change the source RA, Dec to match the center of the observation
    SourceRA,r,a,45.0,0,360,"Source RA (degrees)"
    SourceDEC,r,a,30.0,-90,90,"source DEC (degrees)"

    # The source type should be "SIMPUT"
    SourceType,s,a,"SIMPUT","POINT|GAUSS|IMAGE|LINE|BETA|RAYFILE|DISK|USER|SAOSAC|SIMPUT",,"source"

    # Pointers to your SIMPUT file and the location of the SIMPUT library
    S-SIMPUT-Source,f,a,"sloshing_events_simput.fits",,,"Filename of SIMPUT Catalog"
    S-SIMPUT-Library,f,a,"/usr/local/simput-2.1.2/lib/libsimput.dylib",,,"Path to dynamically linked file libsimput.so"

    # Pointing RA and Dec is up to you, but should be near the source
    RA_Nom,r,a,45.,,,"RA_NOM for dither (degrees)"
    Dec_Nom,r,a,30.,,,"DEC_NOM for dither (degrees)"
    Roll_Nom,r,a,0.,,,"ROLL_NOM for dither (degrees)"

SIMX
++++

Here is an example set of SIMX commands that uses a SIMPUT file made with
pyXSIM:

.. code-block:: bash

    #!/bin/bash
    heainit
    simxinit
    
    punlearn simx
    pset simx mode=hl
    pset simx Exposure=1.0e4
    pset simx UseSimput=yes
    pset simx MissionName=XraySurveyor InstrumentName=HDXI
    pset simx ScaleBkgnd=0.0
    pset simx RandomSeed=24
    
    pset simx SimputFile=spiral_242959_noshift_xrs_simput.fits
    pset simx PointingRA=30.0 PointingDec=45.0
    pset simx OutputFileName=spiral_242959_noshift_xrs
    simx

SOXS
++++

Here is an example set of SOXS commands that uses a SIMPUT file made with
pyXSIM:

.. code-block:: python

    from soxs import instrument_simulator
    simput_file = "snr_simput.fits" # SIMPUT file to be read
    out_file = "evt_mucal.fits" # event file to be written
    exp_time = 30000. # The exposure time in seconds
    instrument = "mucal" # short name for instrument to be used
    sky_center = [30., 45.] # RA, Dec of pointing in degrees
    instrument_simulator(simput_file, out_file, exp_time, instrument,
                         sky_center, overwrite=True)

Refer to the relevant documentation for all of those packages for more details, as well 
as the :ref:`simput` section of the :class:`~pyxsim.event_list.EventList` documentation.