.. _thermal-sources:

Thermal Sources
---------------

X-ray emission from theraml sources in pyXSIM can be generated from models 
which assume collisional ionization is dominant (whether in equilibrium or not)
or a combination of collisional and photoionization processes are relevant. In
the former case, the emission is a function of temperature :math:`T` and 
metallicity :math:`Z`, and is proportional to the density squared:

.. math::

    \varepsilon(E) = n_en_H\Lambda(T, Z, E)

In the case where photoionization is important, the emission is directly 
dependent on the number density of hydrogen:

.. math::

    \varepsilon(E) = n_en_H\Lambda(nH, T, Z, E)

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

By default, setting up a :class:`~pyxsim.source_models.ThermalSourceModel` 
object requires the following arguments:

* ``spectral_model``: The thermal spectral model to assume. Can be a string or 
  :class:`~pyxsim.spectral_models.SpectralModel` instance. Currently, the only
  string value built into pyXSIM is ``"apec"``. 
* ``emin``: The minimum energy for the spectrum in keV.
* ``emax``: The maximum energy for the spectrum in keV.
* ``nchan``: The number of channels in the spectrum. If one is thermally 
  broadening lines (the default), it is recommended that this number create an 
  energy resolution per channel of roughly 1 eV.
* ``Zmet``: The metallicity. Either a floating-point number for a constant
  metallicity, or the name of a yt field for a spatially-varying metallicity.

So creating a default instance is rather simple:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel("apec", 0.1, 11.0, 10000, 0.3)

However, this model is very customizable. There are a number of other optional 
parameters which can be set:

* ``temperature_field``: The yt field to use as the temperature. Must have 
  dimensions of temperature. The default is ``("gas", "temperature")`` for 
  grid-based datasets and ``("PartType0", "Temperature")`` or 
  ``("io", "temperature")`` for particle-based datasets, depending on which is
  available.
* ``emission_measure_field``: The yt field to use as the emission measure. Must
  have dimensions of number density or per-volume. The default is 
  ``("gas", "emission_measure")`` for grid-based datasets. For particle-based 
  datasets, a new field is constructed, using the default density and mass 
  fields of the dataset, and the fields ``("PartType0", "ElectronAbundance")``
  ``("PartType0", "NeutralHydrogenAbundance")`` to construct the electron and
  hydrogen ion number densities if they are present in the dataset.
* ``kT_min``: The minimum temperature in units of keV in the set of temperature
  bins. Default is 0.025.
* ``kT_max``: The maximum temperature in units of keV in the set of temperature
  bins. Default is 64.0.
* ``n_kT``: The number of temperature bins to use. Default is 10000.
* ``kT_scale``: The scaling of the temperature bins, either "linear" or "log".
  Default: "linear"
* ``method``: The method used to generate the photon energies from the spectrum.
  Either ``"invert_cdf"``,
  which inverts the cumulative distribution function of the spectrum, or 
  ``"accept_reject"``, which uses the acceptance-rejection method on the 
  spectrum. The first method should be sufficient for most cases.
* ``thermal_broad``: A boolean specifying whether or not the spectral lines
  should be thermally broadened. Default: True
* ``model_root``: A path specifying where the model files are stored. If not 
  provided, a default location known to pyXSIM is used.
* ``model_vers``: The version identifier string for the model files, e.g. 
  "2.0.2". The default depends on the model used.
* ``var_elem``: Used to specify the abundances of specific elements, whether via
  floating-point numbers of yt fields. A dictionary of elements and values 
  should be specified. See :ref:`var-abund` below for more details.
* ``nolines``: If set to ``True``, the photons for this source will be generated 
  assuming no emission lines. Default: ``False``
* ``abund_table``: The solar abundance table assumed for the different elements.
  See the discussion in :ref:`solar-abund-tables` below for more details. 
  Default: ``"angr"``
* ``prng``: A pseudo-random number generator. Typically will only be specified
  if you have a reason to generate the same set of random numbers, such as for a 
  test or a comparison. Default is the :mod:`numpy.random` module, but a 
  :class:`~numpy.random.RandomState` object or an integer seed can also be used. 

Tweaking the Temperature Bins
+++++++++++++++++++++++++++++

As mentioned above, :class:`~pyxsim.source_models.ThermalSourceModel` bins the 
dataset's cells/particles into a 1-D table of temperatures, each bin containing
a spectrum. It is important that this temperature binning faithfully reflects 
the temperature distribution within the dataset adequately. It may be necessary
to tweak the number, limits, or scaling of the temperature bins. Some example 
situations where it may be necessary to do this are:

* A situation in which there is a lot of low-temperature, high-density gas that 
  is not expected to emit X-rays, in which case one could set ``kT_min`` to a 
  higher value than these temperatures. 
* A situation in which the temperatures in the dataset span a small dynamic 
  range, in which case one would set both ``kT_min`` and ``kT_max`` to bracket 
  this range, and set ``n_kT`` to ensure that the bins are finely spaced. 
* A situation with both low and high temperature gas which are expected to emit 
  X-rays, requiring resolution over a large dynamic range. One could set 
  ``n_kT`` to a large value, or alternatively one could set ``kT_scale="log"`` 
  to adopt logarithmic binning. 

Some degree of trial and error may be necessary to determine the correct setup 
of the temperature bins.

.. _solar-abund-tables:

Changing the Solar Abundance Table
++++++++++++++++++++++++++++++++++

The abundance parameters discussed so far assume abundance of a particular 
element or a number of elements relative to the Solar value. Underlying this
are the values of the Solar abundances themselves. It is possible to change the
Solar abundance table in pyXSIM via the optional ``abund_table`` argument to 
:class:`~pyxsim.source_models.ThermalSourceModel`. By default, pyXSIM assumes 
the `Anders & Grevesse 1989 <http://adsabs.harvard.edu/abs/1989GeCoA..53..197A>`_ 
abundances corresponding to a setting of ``"angr"`` for this parameter, but it 
is possible to use other tables of solar abundances. The other tables included 
which can be used are:

* ``"aspl"``: `Asplund et al. 2009 <http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A>`_
* ``"wilm"``: `Wilms et al. 2000 <http://adsabs.harvard.edu/abs/2000ApJ...542..914W>`_
* ``"lodd"``: `Lodders 2003 <http://adsabs.harvard.edu/abs/2003ApJ...591.1220L>`_

The Solar abundance table can be changed like this:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel("apec", 0.1, 20.0, 10000, 
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

    thermal_model = pyxsim.ThermalSourceModel("apec", 0.1, 20.0, 10000, 
                                              prng=25, abund_table=my_abund)

.. _var-abund:

Variable Abundances
+++++++++++++++++++

By default, :class:`~pyxsim.source_models.ThermalSourceModel` assumes all 
abundances besides H, He, and the trace elements are set by the single value or
yt field provided by the ``Zmet`` parameter. However, more fine-grained control
is possible. :class:`~pyxsim.source_models.ThermalSourceModel` accepts a 
``var_elem`` optional argument to specify which elements should be allowed to
vary freely. The syntax is the same as for ``Zmet``, in that each element set 
can be a single floating-point value or a yt field name corresponding to a field
in the dataset. ``var_elem`` should be a dictionary of key, value pairs where 
the key is the standard abbreviation for the element and the value is the single 
number or field name:

.. code-block:: python

    # Setting abundances by yt field names
    Zmet = ("gas", "metallicity")
    var_elem = {"O": ("gas", "O_fraction"), "Ca": ("gas","Ca_fraction")} 
    source_model = pyxsim.CIESourceModel(0.05, 50.0, 10000, Zmet, var_elem=var_elem)
    
.. code-block:: python

    # Setting abundances by numbers
    Zmet = 0.3
    var_elem = {"O": 0.4, "Ca": 0.5} 
    source_model = pyxsim.ThermalSourceModel(0.05, 50.0, 10000, Zmet, var_elem=var_elem)

Whatever elements are not specified here are assumed to be set as normal, 
whether they are H, He, trace elements, or metals covered by the ``Zmet`` 
parameter. 

.. _nei:

Non-Equilibrium Ionization
++++++++++++++++++++++++++

pyXSIM 2.2.0 and afterward has support for non-equilibrium ionization (NEI) 
emitting plasmas in :class:`~pyxsim.source_models.NEISourceModel`. First, 
one must create a dictionary mapping elements in their different ionization 
states to the corresponding fields in your dataset as seen from yt:

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

    source_model = pyxsim.NEISourceModel(0.3, 1.7, 1000, var_elem=var_elem)

Note that no other elements will be modeled except those which are specified
in ``var_elem``.

Examples
++++++++

Here, we will show several examples of constructing 
:class:`~pyxsim.source_models.ThermalSourceModel` objects. 

An example where we use the default parameters, and a constant 
metallicity:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel("apec", 0.1, 20.0, 10000, 0.5)

An example where we use a metallicity field and change the temperature field:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel("apec", 0.1, 20.0, 10000, 
                                              ("gas", "metallicity"),
                                              temperature_field=("hot_gas","temperature")

An example where we change the limits and number of the temperature bins:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel("apec", 0.1, 20.0, 10000, 0.3,
                                              kT_min=0.1, kT_max=100.,
                                              n_kT=50000)
                                              
An example where we turn off thermal broadening of spectral lines, specify a
directory to find the model files, and specify the model version:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel("apec", 0.1, 20.0, 10000, 0.3,
                                              thermal_broad=False, 
                                              model_root="/Users/jzuhone/data",
                                              model_vers="3.0.3")

An example where we specify a random number generator:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel("apec", 0.1, 20.0, 10000, 0.3,
                                              prng=25)

Turning off line emission:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel("apec", 0.1, 20.0, 10000, 0.3,
                                              prng=25, nolines=True)
