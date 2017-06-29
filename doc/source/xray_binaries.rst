.. _xray-binaries:

Generating Photons from X-ray Binaries
======================================

pyXSIM provides a pair of functions for generating simulated X-ray photons
which originate from X-ray binaries (XRBs) within galaxies. These functions take
simulation datasets with star particles possessing age and metallicity fields
and generates XRB particles via a semi-analytic prescription. Low-mass XRBs (LMXBs)
and high-mass XRBs (HMXBs) are handled seperately in this framework. This is
an experimental feature which is currently in beta.

How it Works
------------

The first step is to take the star particles and determine what the luminosity of
LMXBs and HMXBs is for each one. This is done using the luminosity distribution 
functions from Figure 2 of `Fragos, T., Lehmer, B., Tremmel, M., et al. 2013, ApJ, 764, 41 <http://adsabs.harvard.edu/abs/2013ApJ...764...41F>`_, 
which are functions of stellar population age and metallicity. To convert the
bolometric luminosities to a reference band of 2-10 keV, we use the bolometric 
correction factors from their Table 2.

Given XRB luminosities for each stellar particle, we now need to generate their
associated populations of LMXBs and HMXBs. We generate the LMXBs using the 
probability distribution function :math:`dN/dL` given by Equation 8
of `Gilfanov, M. 2004, MNRAS, 349, 146 <http://adsabs.harvard.edu/abs/2004MNRAS.349..146G>`_,
and for the HMXBs we use Equation 18 of
`Mineo, S., Gilfanov, M., & Sunyaev, R. 2012, MNRAS, 419, 2095 <http://adsabs.harvard.edu/abs/2012MNRAS.419.2095M>`_.

Each XRB particle generated inherits the position and velocity of its host stellar
particle, but is given a small random offset in position from the stellar particle
position. 

Finally, given new a dataset of the XRB particles, we can use this dataset to 
generate a :class:`~pyxsim.photon_list.PhotonList` assuming power-law emission
functions for each XRB particle. LMXBs are given a photon index of :math:`\alpha = 1.56`,
and HMXBs are given a photon index of :math:`\alpha = 2`.

Generating the XRB Particle Dataset
-----------------------------------

The first function, :func:`~pyxsim.source_generators.xray_binaries.make_xrb_particles`,
takes a dataset with star particles and generates a new dataset (understandable by 
yt, and therefore pyXSIM) of XRB particles. We need to create a data object
such as a sphere or rectangular region, and supply to
:func:`~pyxsim.source_generators.xray_binaries.make_xrb_particles` this object
as well as fields for the star particle metallicity and age. We also need
to provide a ``scale_length`` parameter, since the positions of the XRB particles
will be inherited from the star particles but we want to spread them around a bit
in the vicinity of the star particle. If the star particles have a smoothing
length field, that is most ideal, but a (value, unit) tuple such as ``(1.0, "kpc")``
or a :class:`~yt.units.yt_array.YTQuantity` may also be supplied. This example uses
a GIZMO FIRE dataset of a galaxy available from http://yt-project.org/data:

.. code-block:: python

    import yt
    import pyxsim
    
    # Load the dataset
    ds = yt.load("FIRE_M12i_ref11/snapshot_600.hdf5")

    # This dataset has a field for the stellar formation time
    # (in terms of the scale factor), but not the age of the
    # particle, so I must create a field for the particle age
    def _age(field, data):
        # Convert to redshift from scale factor
        z_s = 1.0 / data["PartType4", "StellarFormationTime"] - 1.0
        # Use the dataset's cosmology to determine the age
        age = data.ds.cosmology.t_from_z(0.0) - data.ds.cosmology.t_from_z(z_s)
        age.convert_to_units("Gyr")
        return age
    ds.add_field(("PartType4", "particle_age"), function=_age, units="Gyr", 
                 particle_type=True)

    # Create a sphere object centered on the maximum density
    sp = ds.sphere("max", (0.25, "Mpc"))

    metallicity_field = ("PartType4", "Metallicity") # Metallicity field
    age_field = ("PartType4", "particle_age") # Age field
    scale_length = (1.0, "kpc") # No smoothing length here, so providing
                                # a scalar
                                
    # Create the dataset
    new_ds = make_xrb_particles(sp, metallicity_field, age_field,
                                scale_length)

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
    redshift = 0.01 # Original dataset had z = 0, so putting it out just a bit
    emin = 0.1 # in keV
    emax = 10.0 # in keV
    photons = make_xrb_photons(ds, area, exp_time, redshift, 
                               emin, emax, center=sp.center, 
                               cosmology=ds.cosmology)

Here, we also used the center of the sphere ``sp`` we created earlier as well
as the :class:`~yt.utilities.cosmology.Cosmology` object from the original dataset.
All of the XRB particles are used in the creation of the photons. 

We can then use this :class:`~pyxsim.photon_list.PhotonList` in the usual
pyXSIM way to create a mock observation. The figure below shows an example 
simulation of XRBs, with projected stellar density on the left and the X-ray
image (including the thermal emission from the hot gas of the galaxy) on the 
right.