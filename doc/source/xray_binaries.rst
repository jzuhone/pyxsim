.. _xray-binaries:

Generating Photons from X-ray Binaries
======================================

Generating the XRB Particle Dataset
-----------------------------------

:func:`~pyxsim.source_generators.xray_binaries.make_xrb_particles`

.. code-block:: python

    metallicity_field = ("PartType0", "Metallicity")
    age_field = ("PartType0", "particle_age")
    new_ds = make_xrb_particles(sp, metallicity_field, age_field,
                                scale_length, output_lums="lums.dat")

This returns a yt "in-memory" particle dataset which functions in the 
same way as any other dataset understandable by pyXSIM. 

Generating the XRB Photons
--------------------------

For the second step, to create a :class:`~pyxsim.photon_list.PhotonList` 
instance from this dataset, the helper function 
:func:`~pyxsim.source_generators.xray_binaries.make_xrb_photons` has
been provided. This function is simply a wrapper for the creation of the 
power-law source model for the XRB particles and the creation of the 
:class:`~pyxsim.photon_list.PhotonList` itself:

.. code-block:: python

    area = (20000.0, "cm**2")
    exp_time = (500.0, "ks")
    redshift = 0.01
    emin = 0.1 # in keV
    emax = 10.0 # in keV
    photons = make_xrb_photons(ds, area, exp_time, redshift, emin, emax)

All of the XRB particles are used in the creation of the photons. 

References
----------

The distribution of LMXB and HMXB luminosities as a function of stellar population
age and metallicity is determined from Figure 2 of 
`Fragos, T., Lehmer, B., Tremmel, M., et al. 2013, ApJ, 764, 41
 <http://adsabs.harvard.edu/abs/2013ApJ...764...41F>`_, and the bolometric corrections
are determined from their Table 2.

The probability distribution function :math:`dN/dL` for the LMXBs is given by Equation 8
of `Gilfanov, M. 2004, MNRAS, 349, 146 <http://adsabs.harvard.edu/abs/2004MNRAS.349..146G>`_,
and that of the HMXBs is given by Equation 18 of
`Mineo, S., Gilfanov, M., & Sunyaev, R. 2012, MNRAS, 419, 2095 <http://adsabs.harvard.edu/abs/2012MNRAS.419.2095M>`_.
