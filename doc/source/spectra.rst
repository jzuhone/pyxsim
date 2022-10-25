.. _xray-spectra:

Generating Spectra from Data Objects
====================================

If you have a yt data object (such as a sphere, box, disk) and a source model 
of any sort, then you can also generate spectra from the entire object. This can 
be done in two modes--either in the rest frame of the source, in which case the 
spectrum will be a "count rate" spectrum in units of :math:`\rm{counts}~\rm{s}^{-1}~\rm{keV}^{-1}`, 
or in an observer frame at some distance in which case the spectrum will be in units of 
:math:`\rm{counts}~\rm{cm}^{-2}~\rm{s}^{-1}~\rm{keV}^{-1}`.

Assuming one has a dataset and (say) a sphere object, you can generate spectra 
like this:

.. code-block:: python

    # Here's a power-law source model, but any source model will do

    emin = 0.5
    emax = 40.0
    e0 = 1.0
    alpha
    norm_field = ("gas", "power_law_emission")

    plaw_model = pyxsim.PowerLawSourceModel(e0, emin, emax, norm_field, alpha)

    # Make a count rate spectrum in the source frame
    
    nbins = 1000    

    spec_src = plaw_model.make_spectrum(sp, emin, emax, nbins)

The resulting spectrum ``spec_src`` is a :class:`~soxs.spectrum.CountRateSpectrum` 
object, which has a number of methods in SOXS that can be used to analyze and visualize
it. 

If we instead want to find a spectrum of a source measured by an observer at a specific
distance, then we can specify either a redshift and a cosmology (which uniquely specifies
a distance) or we can set a distance explicitly for a local source.

Here's an example where only a redshift is specified. In this case, a cosmology is assumed 
by default from yt, usually the one associated with the dataset:

.. code-block:: python

    # Make a flux spectrum in the observer frame at some redshift
    
    emin_obs = 2.0
    emax_obs = 20.0
    redshift = 0.1
    
    spec_obs = plaw_model.make_spectrum(sp, emin, emax, nbins, redshift=redshift)

The resulting spectrum ``spec_obs`` is a :class:`~soxs.spectrum.Spectrum` object, which 
has a number of methods in SOXS that can be used to analyze and visualize it.

If you want to choose a different cosmology, specify a yt 
:class:`~yt.utilities.cosmology.Cosmology` object:

.. code-block:: python

    # Make a flux spectrum in the observer frame at some redshift
    # at a specified cosmology
    
    from yt.utilities.cosmology import Cosmology
    
    cosmo = Cosmology(hubble_constant=0.704, omega_matter=1.0-0.728)

    emin_obs = 2.0
    emax_obs = 20.0
    redshift = 0.1
    
    spec_obs = plaw_model.make_spectrum(sp, emin, emax, nbins, redshift=redshift,
                                        cosmology=cosmo)

You can also simply specify a distance in the ``dist`` keyword argument, if the 
source is local (but note that in this case you cannot specify a redshift at the
same time):

.. code-block:: python

    # Make a flux spectrum in the observer frame at some local distance
    
    emin_obs = 2.0
    emax_obs = 20.0
    
    spec_obs = plaw_model.make_spectrum(sp, emin, emax, nbins, dist=(8.0, "kpc"))
   
.. note::

    At this time, Doppler-shifting of photon energies by motions of the emitting 
    material is not available for the creation of spectra in this mode, but it will 
    be available in a future release.

