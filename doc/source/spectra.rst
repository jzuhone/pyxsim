.. _xray-spectra:

Generating Spectra from Data Objects
====================================

If you have a yt data object (such as a sphere, box, disk) and a source model
of any sort, then you can also generate spectra from the entire object. This can
be done in two modes--either in the rest frame of the source, in which case the
spectrum will be a "count rate" spectrum in units of counts s\ :sup:`-1` keV\ :sup:`-1`,
or in an observer frame at some distance in which case the spectrum will be in
units of counts s\ :sup:`-1` cm\ :sup:`-2` keV\ :sup:`-1`


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
distance in their own reference frame, then we can specify either a redshift and a
cosmology (which uniquely specifies a distance) or we can set a distance explicitly
for a local source.

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

Finally, it is also possible to simulate the Doppler shifting of the spectrum from the
velocity field of the source. This is done by specifying the ``normal`` keyword argument,
which can be either be a string corresponding to one of the coordinate axes
(e.g. "x", "y", or "z"), or a 3-vector of the form ``[a, b, c]`` that specifies the
normal direction to the source:

.. code-block:: python

    emin_obs = 2.0
    emax_obs = 20.0

    # Make a flux spectrum in the observer frame at some local distance with
    # Doppler shifting along the z-axis

    spec_obs1 = plaw_model.make_spectrum(sp, emin, emax, nbins, dist=(8.0, "kpc"),
                                         normal="z")

    # Make a flux spectrum in the observer frame at some redshift with
    # Doppler shifting along an arbitrary direction

    spec_obs2 = plaw_model.make_spectrum(sp, emin, emax, nbins, redshift=0.05,
                                         normal=[3.0, -0.2, 1.0])


.. note::

    Doppler shifting via the ``normal`` keyword is only possible if a distance or redshift
    is specified, since otherwise the spectrum is by definition computed in the rest frame
    of the source.
