.. _instruments:

Producing Realistic Observations Using External Packages
========================================================

If you want to produce a realistic simulation of a particular instrumental 
configuration, pyXSIM provides options for exporting its event lists to 
external packages. For `MARX <http://space.mit.edu/ASC/MARX/>`_ and 
`SIMX <http://hea-www.cfa.harvard.edu/simx/>`_, or for more fine-tuned 
use of `SOXS <http://hea-www.cfa.harvard.edu/~jzuhone/soxs>`_, one can use 
SIMPUT files. 

.. note::

    The old built-in instrument models have been deprecated in favor of 
    exporting to one of these packages.

MARX
----

The MARX version needs to be at least 5.3.1. To use SIMPUT with MARX, one only 
needs to change the following lines in the ``marx.par`` file:

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
----

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
----

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

Refer to the relevant documentation for all of those packages for more details,
as well as the :ref:`simput` section of the :class:`~pyxsim.event_list.EventList`
documentation.