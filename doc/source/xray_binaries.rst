.. _xray-binaries:

Generating Photons from X-ray Binaries
======================================

pyXSIM provides a pair of functions for generating simulated X-ray photons
which originate from X-ray binaries (XRBs) within galaxies. These functions take
simulation datasets with star particles possessing age and mass fields
and generates XRB particles via a semi-analytic prescription. Low-mass XRBs (LMXBs)
and high-mass XRBs (HMXBs) are handled seperately in this framework. This is
an experimental feature which is currently in beta.

How it Works
------------

The first step is to take the star particles and determine the total 
number of LMXBs and HMXBs for the distribution. We generate the LMXBs 
using the probability distribution function :math:`dN/dL` given by
Equation 8 of `Gilfanov, M. 2004, MNRAS, 349, 146 <http://adsabs.harvard.edu/abs/2004MNRAS.349..146G>`_,
and for the HMXBs we use Equation 18 of
`Mineo, S., Gilfanov, M., & Sunyaev, R. 2012, MNRAS, 419, 2095 <http://adsabs.harvard.edu/abs/2012MNRAS.419.2095M>`_.
For the HMXBs, we must take into account the star formation rate, which
we estimate by taking the average mass of stars created within a short
period of time previous to the current simulation time. We then determine
via Poisson sampling how many XRBs each star particle will have given its 
mass (most will have either 0 or 1). For each particle, we then use 
:math:`dN/dL` for each type of XRB to determine their luminosities.

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
as well the field for the star particle age. We also need to provide a 
``scale_length`` parameter, since the positions of the XRB particles will be 
inherited from the star particles but we want to spread them around a bit
in the vicinity of the star particle. If the star particles have a smoothing
length field, that is most ideal, but a (value, unit) tuple such as ``(1.0, "kpc")``
or a :class:`~yt.units.yt_array.YTQuantity` may also be supplied. Finally, to make
an estimate of the star formation rate, the total mass of stars created over a short
period of time previous to the current simulation time is computed. This time is
1.0 Gyr by default, but may be set by suppling a new value to the ``sfr_time_range``
parameter.

This example uses a GIZMO FIRE dataset of a galaxy available from 
http://yt-project.org/data:

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

    age_field = ("PartType4", "particle_age") # Age field
    scale_length = (1.0, "kpc") # No smoothing length here, so providing
                                # a scalar
    sfr_time_range = (2.0, "Gyr") # Recent duration over which to compute
                                  # the star formation rate

    # Create a sphere object centered on the maximum density
    sp = ds.sphere("max", (0.25, "Mpc"))

    # Create the dataset
    new_ds = make_xrb_particles(sp, age_field, scale_length)

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
    photons_xrb = make_xrb_photons(new_ds, area, exp_time, redshift, 
                                   emin, emax, center=sp.center, 
                                   cosmology=ds.cosmology)

    # We now show how to compute photons from the gas, generate events,
    # and add the events together
    
    source_model = pyxsim.ThermalSourceModel("apec", 0.1, 10.0, 10000)

    photons_gas = pyxsim.PhotonList.from_data_source(sp,  redshift, area, exp_time, source_model)

    events_xrb = photons_xrb.project_photons("z", [30.0, 45.0], absorb_model="tbabs", nH=0.02)
    events_gas = photons_gas.project_photons("z", [30.0, 45.0], absorb_model="tbabs", nH=0.02)

    events = events_xrb + events_gas

Here, we also used the center of the sphere ``sp`` we created earlier as well
as the :class:`~yt.utilities.cosmology.Cosmology` object from the original dataset.
All of the XRB particles are used in the creation of the photons. 

We can then use this :class:`~pyxsim.photon_list.PhotonList` in the usual
pyXSIM way to create a mock observation. The figure below shows an example 
simulation of XRBs, with projected stellar density on the left and the 
X-ray image (including the thermal emission from any hot gas in the galaxy)
on the right, produced using the ACIS-I simulator built into SOXS.

.. image:: _images/gizmo_xrbs.prng

Other Examples
--------------

The following two examples use galaxy datasets which can be downloaded 
from http://yt-project.org/data. Only the code necessary to make the 
XRB particles is shown here, the rest is the same. 

An ART dataset:

.. code-block:: python

    import yt
    import pyxsim

    ds = yt.load("sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art")
    
    scale_length = (1.0, "kpc")
    age_field = ("STAR", "age")
    mass_field = ("STAR", "particle_mass")
    
    sp = ds.sphere("max", (0.25, "Mpc"))

    new_ds = make_xrb_particles(sp, age_field, scale_length)

.. image:: _images/art_xrbs.prng

An Enzo dataset:

.. code-block:: python

    import yt
    import pyxsim

    @yt.particle_filter(requires=["particle_type"], filtered_type='all')
    def stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 2
        return filter
    
    ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    ds.add_particle_filter('stars')
    
    scale_length = (1.0, "kpc")
    age_field = ("stars", "age")
    mass_field = ("stars", "particle_mass")

    sp = ds.sphere("max", (0.25, "Mpc"))

    new_ds = make_xrb_particles(sp, age_field, scale_length)

.. image:: _images/enzo_xrbs.prng
