.. _instruments:

Convolving Events with Instrumental Responses
=============================================

Built-In Instrument Simulators
------------------------------

pyXSIM provides the ability to perform "quick-and-dirty" convolutions with approximate
representations of real X-ray instruments. The accuracy of these representations is 
limited, since they assume the following simplifications:

* A square field of view without chip gaps
* The spectral response, PSF, and effective area are position-independent.
* The PSF is assumed to have a Gaussian shape
* No instrumental background is added

If you only need an approximate representation of what an X-ray observation of your source
would look like, using a built-in instrument simulator should be sufficient for your purposes. 

.. warning::

    If you want to export simulated events to an external instrument simulator
    such as MARX or SIMX (see below), DO NOT use any of these built-in instrument 
    simulators within pyXSIM. This will be the job of the instrument simulator itself to do.

Using a Built-In Instrument Simulator
+++++++++++++++++++++++++++++++++++++

pyXSIM comes with the following built-in instrument simulators:

* ``ACIS_I``: ACIS-I Cycle 18 on-axis ARF and RMF, with 0.492" central pixel and 0.5" FWHM PSF
* ``ACIS_S``: ACIS-S Cycle 18 on-axis ARF and RMF, with 0.492" central pixel and 0.5" FWHM PSF
* ``Hitomi_SXS``: Hitomi Soft X-ray Spectrometer ARF and RMF, with 30.64" central pixel and 1.2' FWHM PSF
* ``Athena_WFI``: Athena Wide-Field Imager ARF and RMF, with 2.23" central pixel and 5.0" FWHM PSF
* ``Athena_XIFU``: Athena X-ray Integral Field Unit ARF and RMF, with 4.56" central pixel and 5.0" FWHM PSF
* ``XRS_Imager``: X-ray Surveyor Imager ARF and RMF, with 0.33" central pixel and 0.5" FWHM PSF
* ``XRS_Calorimeter``: X-ray Surveyor Calorimeter ARF and RMF, with 1.03" central pixel and 0.5" FWHM PSF

When an instrument simulator is called, the following operations are applied to the input events, in
this order.

1. The event positions are re-binned from the original simulation pixelization to the one appropriate
   for the detector simulation.
2. The event positions are smoothed using a Gaussian PSF. 
3. Using the effective area curve from the selected ARF, events are selected or rejected for observation.
4. The observed event energies are convolved with the selected RMF to produce the observed energy channels. 

Assuming one has an :class:`~pyxsim.event_list.EventList` object handy, to generate a new event list
passed through one of these :class:`~pyxsim.instruments.InstrumentSimulator` classes, one only needs to call
the instrument simulator with the events as an argument:

.. code-block:: python

    >>> from pyxsim import ACIS_S, Hitomi_SXS
    >>> aciss_events = ACIS_S(events)
    >>> hitomi_events = Hitomi_SXS(events)

The specific effects of the instrument simulator can be turned on or off, which are handled with the
following boolean keyword arguments, all of which default to ``True``:

* ``rebin``: Controls whether or not the events are rebinned.
* ``convolve_psf``: Controls whether or not the PSF smoothing is applied.
* ``convolve_arf``: Controls whether or not the events are selected according to the effective area curve.
* ``convolve_rmf``: Controls whether or not the event energies are convolved with the spectral response. Note that
  if ``convolve_arf=False``, the energies will not be convolved either. 

For example, if one only wanted to detect events using the ARF and convolve the energies with the RMF, one
could do this:

.. code-block:: python

    >>> from pyxsim import ACIS_S
    >>> new_events = ACIS_S(events, rebin=False, convolve_psf=False)

Designing Your Own Instrument Simulator
+++++++++++++++++++++++++++++++++++++++

If you want to design an instrument simulator yourself for use with pyXSIM, it is fairly simple.
You need to provide the following parameters in the call to :class:`~pyxsim.instruments.InstrumentSimulator`, 
in this order: 

* ``dtheta``: The width of the reference (central) pixel in degrees.
* ``nx``: The number of resolution elements (pixels) on a side across the field of view.
* ``psf_scale``: The FWHM of the Gaussian PSF in degrees. 
* ``arf``: The path to the ARF file you want to use. 
* ``rmf``: The path to the RMF file you want to use. 

.. warning::

    It goes without saying that the ARF and RMF you choose MUST be consistent with each other; e.g., 
    have the same energy binning, etc.
    
Here we show an example of how to set up an instrument simulator, using the parameters that pyXSIM
uses for the built-in ACIS-S simulator:

.. code-block:: python

    from pyxsim import InstrumentSimulator

    ACIS_S = InstrumentSimulator(0.0001366667, 8192, 0.0001388889,
                                 "aciss_aimpt_cy18.arf",
                                 "aciss_aimpt_cy18.rmf")

Producing More Realistic Observations Using External Packages
-------------------------------------------------------------

If you want to produce a more realistic simulation of a particular instrumental configuration,
pyXSIM provides options for exporting its event lists to external packages. For 
`MARX <http://space.mit.edu/ASC/MARX/>`_ and `SIMX <http://hea-www.cfa.harvard.edu/simx/>`_,
one can use SIMPUT files. 

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

Refer to the relevant documentation for both of those packages for
more details, as well as the :ref:`simput` section of the :class:`~pyxsim.event_list.EventList`
documentation.