.. _thermal-sources:

Thermal Sources
---------------

X-ray emission from thermal sources in pyXSIM can be generated from models 
which assume collisional ionization is dominant (whether in equilibrium or not)
or a combination of collisional and photoionization processes are relevant. In
the former case, the emission is a function of temperature :math:`T` and 
metallicity :math:`Z`, and is proportional to the density squared:

.. math::

    \varepsilon(E) = n_en_H\Lambda(T, Z; E)

In the case where photoionization is important, the emission is directly 
dependent on the number density of hydrogen:

.. math::

    \varepsilon(E) = n_en_H\Lambda(n_H, T, Z; E)

In either case, the metallicity :math:`Z` may be a single value encompassing 
all metals or a vector corresponding to a list of individual elements. Since 
generating spectra and/or fluxes from each particle or cell would be 
prohibitively expensive, the emission is first binned into a 1-D or 2-D table
(depending on whether or not the spectrum itself is density-dependent or only
dependent on the temperature), and for each bin a spectrum or flux is calculated. 
For each cell or particle, the spectrum or flux is then linearly interpolated 
from this table. The accuracy of this method is sufficient for most purposes. 

.. warning::

    Thermal source models at minimum require a temperature field and a density
    field. In addition, some information regarding the number densities of 
    hydrogen and free electrons need to be provided, so that a field for the
    emission measure :math:`n_en_HdV` can be constructed. This requires the
    ``("gas","H_nuclei_density")`` (total number density of all species of
    hydrogen) and ``("gas","El_number_density")`` (number density of free
    electrons) fields to be present in your dataset, which should be the case
    if your simulation keeps track of individual species. If they are not, 
    and you only have a ``("gas","density")`` field, you may assume full 
    ionization by loading your dataset like this: 
    ``ds = yt.load(filename, default_species_fields="ionized")``. 

Three types of thermal source models are available, which will be described
in turn below. 

.. _cie-source-model:

Thermal Sources in Collisional Ionization Equilbrium (CIE)
==========================================================

The :class:`~pyxsim.source_models.thermal_sources.CIESourceModel` class
simulates thermal emission under the assumption of collisional ionization
equilibrium (CIE). By default, setting up a 
:class:`~pyxsim.source_models.thermal_sources.CIESourceModel` object requires 
the following arguments:

* ``model``: The specific CIE model to use. Options are ``"apec"``, ``"spex"``,
  ``"mekal"``, or ``"cloudy"``. 
* ``emin``: The minimum energy for the spectrum in keV.
* ``emax``: The maximum energy for the spectrum in keV.
* ``nbins``: The number of bins in the spectrum. 
* ``Zmet``: The metallicity. Either a floating-point number for a constant
  metallicity, or the name of a yt field for a spatially-varying metallicity.

Exactly what abundances are set by the ``Zmet`` parameter depends on the 
model used:

* ``"apec"`` and ``"spex"``: Includes C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, 
  Fe, Ni (He fixed at abundance table value, Li, Be, B, F, Na, P, Cl, K, 
  Sc, Ti, V, Cr, Mn, Co, Cu, Zn fixed at solar).
* ``"mekal"``: Includes C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni 
  (He fixed at abundance table value, other elements not modeled)
* ``"cloudy"``: He fixed at abundance table value, all higher elements up 
  Zn to included.

Creating a default instance is rather simple:

.. code-block:: python

    thermal_model = pyxsim.CIESourceModel("apec", 0.1, 11.0, 10000, 0.3)

However, this model is very customizable. There are a number of other optional 
parameters which can be set:

* ``binscale``: The scale of the energy binning: "linear" or "log". 
  Default: "linear"
* ``temperature_field``: The yt field to use as the temperature. Must have 
  dimensions of temperature. The default is ``("gas", "temperature")``.
* ``emission_measure_field``: The yt field to use as the emission measure. Must
  have dimensions of number density or per-volume. The default is 
  ``("gas", "emission_measure")``. 
* ``h_fraction``: The hydrogen mass fraction. If a float, assumes a constant 
  mass fraction of hydrogen throughout. If a string or tuple of strings, 
  is taken to be the name of the hydrogen fraction field from yt. Defaults to
  the appropriate value for the chosen abundance table.
* ``kT_min``: The minimum temperature in units of keV. Default is 0.025.
* ``kT_max``: The maximum temperature in units of keV. Default is 64.0.
* ``max_density``: The maximum mass density of the cells or particles to use 
  when generating photons. If a float, the units are assumed to be g/cm**3. 
  can also be a ``YTQuantity`` or a float, string tuple such as 
  ``(5.0e-25, "g/cm**3")`` Default: None, meaning no maximum density.
* ``method``: The method used to generate the photon energies from the spectrum.
  Either ``"invert_cdf"``,
  which inverts the cumulative distribution function of the spectrum, or 
  ``"accept_reject"``, which uses the acceptance-rejection method on the 
  spectrum. The first method should be sufficient for most cases.
* ``thermal_broad``: A boolean specifying whether or not the spectral lines
  should be thermally broadened. Only available for the ``"apec"`` and 
  ``"spex"`` models. Default: True
* ``model_root``: A path specifying where the model files are stored. If not 
  provided, a default location known to pyXSIM is used.
* ``model_vers``: The version identifier string for the model files, e.g. 
  "2.0.2". The default depends on the model used. Currently only implemented
  for the ``"apec"`` or ``"spex"`` models.
* ``var_elem``: Optionally used to specify the abundances of specific elements, 
  whether via floating-point numbers or yt fields. A dictionary of elements and 
  values should be specified. See :ref:`var-abund` below for more details.
* ``nolines``: If set to ``True``, the photons for this source will be generated 
  assuming no emission lines. Only available for the ``"apec"`` and ``"spex"`` 
  models. Default: ``False``
* ``abund_table``: The solar abundance table assumed for the different elements.
  See the discussion in :ref:`solar-abund-tables` below for more details. 
  Default: ``"angr"``
* ``prng``: A pseudo-random number generator. Typically will only be specified
  if you have a reason to generate the same set of random numbers, such as for a 
  test or a comparison. Default is the :mod:`numpy.random` module, but a 
  :class:`~numpy.random.RandomState` object or an integer seed can also be used. 

.. _solar-abund-tables:

Changing the Solar Abundance Table
++++++++++++++++++++++++++++++++++

The abundance parameters discussed so far assume abundance of a particular 
element or a number of elements relative to the Solar value. Underlying this
are the values of the Solar abundances themselves. It is possible to change the
Solar abundance table in pyXSIM via the optional ``abund_table`` argument to 
:class:`~pyxsim.source_models.thermal_sources.CIESourceModel`. By default, 
pyXSIM assumes the `Anders & Grevesse 1989 <http://adsabs.harvard.edu/abs/1989GeCoA..53..197A>`_ 
abundances corresponding to a setting of ``"angr"`` for this parameter, but it 
is possible to use other tables of solar abundances. tables included 
which can be used are:

* ``"angr"``: `Anders & Grevesse 1989 <http://adsabs.harvard.edu/abs/1989GeCoA..53..197A>`_
* ``"aspl"``: `Asplund et al. 2009 <http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A>`_
* ``"wilm"``: `Wilms et al. 2000 <http://adsabs.harvard.edu/abs/2000ApJ...542..914W>`_
* ``"lodd"``: `Lodders 2003 <http://adsabs.harvard.edu/abs/2003ApJ...591.1220L>`_
* ``"feld"``: `Feldman 1992 <https://ui.adsabs.harvard.edu/abs/1992PhyS...46..202F>`_
* ``"cl17.03"``: The abundances used by default in Cloudy 17.03.

The Solar abundance table can be changed like this:

.. code-block:: python

    thermal_model = pyxsim.CIESourceModel("apec", 0.1, 20.0, 10000, 
                                          ("gas","metallicity"),
                                          prng=25, abund_table='lodd')

Alternatively, one can supply their own abundance table by providing a NumPy 
array, list, or tuple of abundances 30 elements in length corresponding to the
Solar abundances relative to hydrogen in the order of H, He, Li, Be, B, C, N, O,
F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, 
and Zn. An example:

.. code-block:: python

    my_abund = np.array([1.00E+00, 8.51E-02, 1.12E-11, 2.40E-11, 5.01E-10,
                         2.69E-04, 6.76E-05, 4.90E-04, 3.63E-08, 8.51E-05,
                         1.74E-06, 3.98E-05, 2.82E-06, 3.24E-05, 2.57E-07,
                         1.32E-05, 3.16E-07, 2.51E-06, 1.07E-07, 2.19E-06,
                         1.41E-09, 8.91E-08, 8.51E-09, 4.37E-07, 2.69E-07,
                         3.16E-05, 9.77E-08, 1.66E-06, 1.55E-08, 3.63E-08])

    thermal_model = pyxsim.CIESourceModel("spex", 0.1, 20.0, 10000, 
                                          prng=25, abund_table=my_abund)

.. note:: 

    Currently the solar abundance table cannot be changed for the ``"cloudy"``
    model. It is set to ``"feld"``. 

.. _var-abund:

Variable Abundances
+++++++++++++++++++

As noted above, by default :class:`~pyxsim.source_models.CIESourceModel` assumes 
all abundances besides H, He, and perhaps some trace elements are set by the single
value or yt field provided by the ``Zmet`` parameter. However, more fine-grained 
control is possible. :class:`~pyxsim.source_models.CIESourceModel` accepts a 
``var_elem`` optional argument to specify which elements should be allowed to vary
freely. The syntax is the same as for ``Zmet``, in that each element set can be a 
single floating-point value or a yt field name corresponding to a field in the 
dataset. ``var_elem`` should be a dictionary of key, value pairs where the key is 
the standard abbreviation for the element and the value is the single number or 
field name:

.. code-block:: python

    # Setting abundances by yt field names
    Zmet = ("gas", "metallicity")
    var_elem = {"O": ("gas", "O_fraction"), "Ca": ("gas","Ca_fraction")} 
    source_model = pyxsim.CIESourceModel("cloudy", 0.05, 50.0, 10000, Zmet, var_elem=var_elem)
    
.. code-block:: python

    # Setting abundances by numbers
    Zmet = 0.3
    var_elem = {"O": 0.4, "Ca": 0.5} 
    source_model = pyxsim.CIESourceModel("mekal", 0.05, 50.0, 10000, Zmet, var_elem=var_elem)

Whatever elements are not specified here are assumed to be set as normal, 
whether they are H, He, trace elements, or metals covered by the ``Zmet`` 
parameter. The abundances that you can specify in ``var_elem`` depend on 
the model being used:

* ``"apec"`` and ``"spex"``: Can vary any element He and higher up to Zn
* ``"mekal"``: Can vary He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni 
* ``"cloudy"``: Can vary C, N, O, Ne, Fe, S, Si, Ca, and Mg

Examples
++++++++

Here, we will show several examples of constructing 
:class:`~pyxsim.source_models.thermal_models.CIESourceModel` objects. 

An example where we use the default parameters, and a constant 
metallicity:

.. code-block:: python

    thermal_model = pyxsim.CIESourceModel("apec", 0.1, 20.0, 10000, 0.5)

An example where we use a metallicity field and change the temperature field:

.. code-block:: python

    thermal_model = pyxsim.CIESourceModel("apec", 0.1, 20.0, 10000, 
                                          ("gas", "metallicity"),
                                          temperature_field=("hot_gas","temperature")

An example where we change the limits of the temperature, and use the MeKaL
model:

.. code-block:: python

    thermal_model = pyxsim.CIESourceModel("mekal", 0.1, 20.0, 10000, 0.3,
                                          kT_min=0.1, kT_max=100.)
                                              
An example where we turn off thermal broadening of spectral lines, specify a
directory to find the model files, and specify the model version:

.. code-block:: python

    thermal_model = pyxsim.CIESourceModel("apec", 0.1, 20.0, 10000, 0.3,
                                          thermal_broad=False, 
                                          model_root="/Users/jzuhone/data",
                                          model_vers="3.0.3")

An example where we specify a random number generator and use the Cloudy
model:

.. code-block:: python

    thermal_model = pyxsim.CIESourceModel("cloudy", 0.1, 20.0, 10000, 0.3,
                                          prng=25)

Turning off line emission for the ``"apec"`` model:

.. code-block:: python
    
    thermal_model = pyxsim.CIESourceModel("apec", 0.1, 20.0, 10000, 0.3,
                                          prng=25, nolines=True)

.. _nei:

Non-Equilibrium Ionization
==========================

pyXSIM has support for emission from plasmas in a non-equilibrium ionization
state in the :class:`~pyxsim.source_models.thermal_sources.NEISourceModel`.
In this case, it is assumed that the NEI calculation for the various ionization
states has been carried out in your simulation code, so that you have fields
available for each element and ionization state that you want to generate
emission from. 

To use :class:`~pyxsim.source_models.thermal_sources.NEISourceModel`, one must 
first create a dictionary mapping elements in their different ionization states 
to the corresponding fields in your dataset as seen from yt, or single 
floating-point values. The ionization states in the keys of this dictionary 
are given in the ``"{elem}^{ion}"`` format, where ``ion=0`` is neutral, 
``ion=1`` is singly ionized, and so on. 

Here is an example from a FLASH dataset:

.. code-block:: python

    # The dict mapping ionization states of different elements to different
    # yt fields
    var_elem = {"H^1": ("flash", "h   "),
                "He^0": ("flash", "he  "),
                "He^1": ("flash", "he1 "),
                "He^2": ("flash", "he2 "),
                "O^0": ("flash", "o   "),
                "O^1": ("flash", "o1  "),
                "O^2": ("flash", "o2  "),
                "O^3": ("flash", "o3  "),
                "O^4": ("flash", "o4  "),
                "O^5": ("flash", "o5  "),
                "O^6": ("flash", "o6  "),
                "O^7": ("flash", "o7  "),
                "O^8": ("flash", "o8  ")
               }

Unlike the :class:`~pyxsim.source_models.thermal_sources.CIESourceModel`, for 
the :class:`~pyxsim.source_models.thermal_sources.NEISourceModel` source all
elements and ionizations must be specified in the ``var_elem`` dictionary,
which is now required. There is no separate ``Zmet`` which can be set. The
required arguments are:

* ``emin``: The minimum energy for the spectrum in keV.
* ``emax``: The maximum energy for the spectrum in keV.
* ``nbins``: The number of bins in the spectrum. 
* ``var_elem``: Used to specify the abundances of specific elements, whether 
  via floating-point numbers or yt fields. A dictionary of elements and values 
  should be specified. 

All other optional keyword arguments are the same as in the 
:class:`~pyxsim.source_models.thermal_sources.CIESourceModel`, see above for
details. The :class:`~pyxsim.source_models.thermal_sources.NEISourceModel`
is currently only compatible with the ``"apec"`` emission model. An example 
invocation is:

.. code-block:: python

    source_model = pyxsim.NEISourceModel(0.3, 1.7, 1000, var_elem)

Note that no other elements will be modeled except those which are specified
in ``var_elem``.

.. _igm-source-model:

IGM Source Model
================

The :class:`~pyxsim.source_models.thermal_sources.IGMSourceModel` is 
a source model for a thermal plasma including photoionization and 
resonant scattering from the CXB, based on 
`Khabibullin & Churazov 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4972K/>`_ 
and `Churazov et al. 2001 <https://ui.adsabs.harvard.edu/abs/2001MNRAS.323...93C/>`_.
Because of the included effects of photoionization and resonant 
scattering, this model is dependent on the hydrogen number density in
an explicit way, aside from the normalization.

This model is appropriate for simulation emission from low-density, 
high-temperature plasmas such as the warm-hot intergalactic medium (WHIM) and
the outskirts of the circumgalactic medium (CGM). The densities and 
temperatures involved are :math:`n_H \sim 10^{-7} - 10^{-4} \rm{cm}^{-3}` and
:math:`T \sim 10^5 - 10^7` K. For resonant scattering, it is assumed that 
a fraction of CXB photons are scattering off of heavy ions, enhancing line
emission. 

For temperatures higher than :math:`kT \sim 1.09` keV, the emission is
essentially density-independent (aside from the normalization) and a 
Cloudy-based CIE model is used to compute the spectrum. This model assumes the
abundance tables from Feldman 1992 (``"feld"``) and currently cannot be changed 
to another.

The arguments for :class:`~pyxsim.source_models.thermal_sources.IGMSourceModel`
are very similar to :class:`~pyxsim.source_models.thermal_sources.CIESourceModel`.
Required arguments are:

* ``emin``: The minimum energy for the spectrum in keV.
* ``emax``: The maximum energy for the spectrum in keV.
* ``nbins``: The number of bins in the spectrum. 
* ``Zmet``: The metallicity. Either a floating-point number for a constant
  metallicity, or the name of a yt field for a spatially-varying metallicity.

For the :class:`~pyxsim.source_models.thermal_sources.IGMSourceModel`, He is 
fixed at abundance table value, and all higher elements up Zn to included in
``Zmet``. Optional arguments are:

* ``binscale``: The scale of the energy binning: "linear" or "log". 
  Default: "linear"
* ``resonant_scattering``: Whether or not to include the effects of resonant 
  scattering from CXB photons. Default: False
* ``cxb_factor``: The fraction of the CXB photons that are resonant scattered 
  to enhance the lines. Default: 0.5
* ``nh_field``: The yt field to use as the number density of hydrogen. 
  Must have number density dimensions. The default is ``("gas", "H_nuclei_density")``.
* ``temperature_field``: The yt field to use as the temperature. Must have 
  dimensions of temperature. The default is ``("gas", "temperature")``.
* ``emission_measure_field``: The yt field to use as the emission measure. Must
  have dimensions of number density or per-volume. The default is 
  ``("gas", "emission_measure")``. 
* ``h_fraction``: The hydrogen mass fraction. If a float, assumes a constant 
  mass fraction of hydrogen throughout. If a string or tuple of strings, 
  is taken to be the name of the hydrogen fraction field from yt. Defaults to
  the appropriate value for the Feldman abundance table.
* ``kT_min``: The minimum temperature in units of keV. Default is 0.00431.
* ``kT_max``: The maximum temperature in units of keV. Default is 64.0.
* ``max_density``: The maximum mass density of the cells or particles to use 
  when generating photons. If a float, the units are assumed to be g/cm**3. 
  can also be a ``YTQuantity`` or a float, string tuple such as 
  ``(5.0e-25, "g/cm**3")`` Default: None, meaning no maximum density.
* ``method``: The method used to generate the photon energies from the spectrum.
  Either ``"invert_cdf"``,
  which inverts the cumulative distribution function of the spectrum, or 
  ``"accept_reject"``, which uses the acceptance-rejection method on the 
  spectrum. The first method should be sufficient for most cases.
* ``var_elem``: Optionally used to specify the abundances of specific elements, 
  whether via floating-point numbers or yt fields. A dictionary of elements and 
  values should be specified. See :ref:`var-abund` below for more details.
* ``prng``: A pseudo-random number generator. Typically will only be specified
  if you have a reason to generate the same set of random numbers, such as for a 
  test or a comparison. Default is the :mod:`numpy.random` module, but a 
  :class:`~numpy.random.RandomState` object or an integer seed can also be used. 

Examples
++++++++

A simple invocation of the IGM model using a single metallicity field, and
log-spaced energy binning: 

.. code-block:: python

    source_model = pyxsim.IGMSourceModel(0.1, 5.0, 1000, 
                                         ("gas","metallicity"), binscale="log")

Turning on resonant scattering, assuming 30% of the CXB photons are scattered:

.. code-block:: python

    source_model = pyxsim.IGMSourceModel(0.1, 5.0, 1000, 
                                         ("gas","metallicity"),
                                         resonant_scattering=True,
                                         cxb_factor=0.3, 
                                         binscale="log")
                                         
Specifying the abundances of C, N, and Fe separately:

.. code-block:: python

    var_elem = {"C": ("gas", "C_fraction"), 
                "N": ("gas", "N_fraction"),
                "Fe": ("gas", "Fe_fraction")}
           
    source_model = pyxsim.IGMSourceModel(0.1, 5.0, 1000, 
                                         ("gas","metallicity"),
                                         resonant_scattering=True,
                                         cxb_factor=0.3, 
                                         binscale="log", 
                                         var_elem=var_elem)

.. _hot-gas-filter:

Filtering Out Non-X-ray Emitting Gas
====================================

In simulations where may gas phases are present, there may be a significant
amount of thermal gas that is not expected to be emitting in X-rays. For any
of the thermal source models detailed above, there are various ways to ensure
that this gas is not operated on in pyXSIM. 

pyXSIM-based Filtering
++++++++++++++++++++++

The first is to make use of the ``kT_min`` and ``kT_max`` keyword arguments:

.. code-block:: python

    thermal_model = pyxsim.CIESourceModel("apec", 0.1, 20.0, 10000, 
                                          ("gas","metallicity"),
                                          kT_min=0.1,
                                          kT_max=50.0,
                                          prng=25, abund_table='lodd')

where both ``kT_min`` and ``kT_max`` are in units of keV. It may also be useful
to specify a maximum density above which no emission should be calculated with the
``max_density`` keyword argument:

.. code-block:: python

        source_model = pyxsim.IGMSourceModel(0.1, 5.0, 1000, 
                                             ("gas","metallicity"), 
                                             max_density=(5.0e-25, "g/cm**3"),
                                             binscale="log")

yt-based Filtering
++++++++++++++++++

If you want more detailed control over which cells or particles may get used,
then you need to use one of 
`various methods in yt for dataset filtering <https://yt-project.org/doc/analyzing/filtering.html>`_. 

AMR cell-based Filtering
^^^^^^^^^^^^^^^^^^^^^^^^

For example, if your dataset is AMR cell-based, then the use of a 
`yt cut region <https://yt-project.org/doc/analyzing/filtering.html#cut-regions>`_ is recommended. 
In this case, we exclude all gas above :math:`T = 3 \times 10^5 \rm{K}`, below 
:math:`\rho = 5 \times 10^{-25}~\rm{g}~\rm{cm}^{-3}`, and include no gas with star formation. 

.. code-block:: python

    # this example takes a box region and makes cuts on density, temperature, and star
    # formation rate
    
    c = ds.find_min(("gas", "gravitational_potential")) # center of box
    w = ds.quan(1.0, "Mpc") # width of box
    le = c-0.5*w # left edge of box
    re = c+0.5*w # right edge of box
    box = ds.box(le, re) # create the box
    
    # chain these conditions together
    hot_box = box.include_above(("gas", "temperature"), 3e5, "K")
    hot_diffuse_box = hot_box.include_below(("gas", "density"), 5e-25), "g/cm**3")
    xray_box = hot_diffuse_box.include_equal(("gas", "star_formation_rate"), 0.0), "Msun/yr")

This is a new data container exactly like a sphere or box object that can be used to 
create photons, make fields, etc.

Particle-based Filtering
^^^^^^^^^^^^^^^^^^^^^^^^

In the case of particle data (including particle-ish data like Arepo Voronoi cells), 
it makes most sense to use a 
`yt particle filter <https://yt-project.org/doc/analyzing/filtering.html#filtering-particle-fields>`_.
This creates a new particle type that can be used in analysis in the same way as
regular particle types. 

Here is an example of instantiating a particle filter for an Arepo dataset.
In this case, we exclude all gas above :math:`T = 3 \times 10^5 \rm{K}`, below 
:math:`\rho = 5 \times 10^{-25}~\rm{g}~\rm{cm}^{-3}`, and include no gas with star 
formation. 

.. code-block:: python

    # define hot gas filter 
    def hot_gas(pfilter, data):
        pfilter1 = data[pfilter.filtered_type, "temperature"] > 3.0e5
        pfilter2 = data["PartType0", "StarFormationRate"] == 0.0
        pfilter3 = data[pfilter.filtered_type, "density"] < 5e-25
        return pfilter1 & pfilter2 & pfilter3
    # add the filter to yt itself
    yt.add_particle_filter("hot_gas", function=hot_gas,
                           filtered_type='gas', requires=["temperature","density"])
    
    # load dataset and assign filter
    ds = yt.load("cutout_136.hdf5")
    ds.add_particle_filter("hot_gas")

Note that for this dataset the ``"gas"`` and ``"PartType0"`` field types are the
same.