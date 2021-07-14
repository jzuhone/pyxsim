.. _source-models:

Source Models for Generating Photons
====================================

pyXSIM comes with three pre-defined ``SourceModel`` types for generating a new
:class:`~pyxsim.photon_list.PhotonList`, for use with the 
:func:`~pyxsim.make_photons` function. Though these should cover the vast 
majority of use cases, there is also the option to design your own source model. 

.. _thermal-sources:

Thermal Sources
---------------

:class:`~pyxsim.source_models.ThermalSourceModel` assumes the emission of a hot 
thermal plasma can be described by a model that is only dependent on temperature 
and metallicity, and is proportional to the density squared:

.. math::

    \varepsilon(E) = n_en_H\Lambda(T, Z, E)

:class:`~pyxsim.source_models.ThermalSourceModel` requires the use of a thermal
spectral model, described in the next sub-section. From this spectral model, 
which depends on temperature and metallicity, a spectrum of photon energies can
be generated from each cell or particle. Since generating a new spectrum for 
each cell would be prohibitively expensive, the cells are binned into a 1-D 
table of temperatures, and for each bin a spectrum is calculated. Provided the
bins are finely spaced enough, the accuracy of this method is sufficient for 
most purposes. 

.. warning::

    This only works if your dataset has a `("gas", "emission_measure")`
    field from yt, which is defined by yt if you have species defined in 
    your dataset such that yt detects them and generates the
    `("gas", "H_nuclei_density")` (total number density of all species of
    hydrogen) and `("gas", "El_number_density")` (number density of free
    electrons) fields. If you do not have these fields defined in your
    dataset, you may assume full ionization by loading your dataset like
    this: ``ds = yt.load(filename, default_species_fields="ionized")``. 

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
    var_elem = {"O": "oxygen", "Ca": "calcium"} 
    source_model = pyxsim.ThermalSourceModel(0.05, 50.0, 10000, Zmet, var_elem=var_elem)
    
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
emitting plasmas in :class:`~pyxsim.source_models.ThermalSourceModel`. First, 
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

Note that no other elements will be modeled except those which are specified
in ``var_elem``.

The flag for NEI must be set ``nei=True`` when making the model object. 
Note that since the NEI tables are not bundled with pyXSIM, they must be 
downloaded from the `AtomDB website <http://www.atomdb.org>`_ and one must
specify their location in ``model_root``. One may also have to change the 
``model_vers`` string if the model version is not the default ``"v3.0.9"``.

.. code-block:: python

    # model files are located here
    model_root = "/Users/jzuhone/atomdb_v3.0.9"

    source_model = pyxsim.ThermalSourceModel("apec", 0.3, 1.7, 1000, 
                                             ("gas","metallicity"), nei=True, 
                                             model_root=model_root,
                                             var_elem=var_elem)


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

.. _power-law-sources:

Power-Law Sources
-----------------

:class:`~pyxsim.source_models.PowerLawSourceModel` assumes that the emission can
be described by a pure power law:

.. math::

    \varepsilon(E) = K\left(\frac{E}{E_0}\right)^{-\alpha}, E_{\rm min} \leq E \leq E_{\rm max}
    
between the energies ``emin`` and ``emax``, with a power-law spectral index 
``alpha``. The power law normalization :math:`K` is represented by an 
``emission_field`` specified by the user, which must have units of counts/s/keV 
in the source rest frame. ``alpha`` may be a single floating-point number 
(implying the spectral index is the same everywhere), or a field specification
corresponding to a spatially varying spectral index. A reference energy ``e0`` 
(see above equation) must also be specified.

Examples
++++++++

An example where the spectral index is the same everywhere:

.. code-block:: python

    e0 = (1.0, "keV") # Reference energy
    emin = (0.01, "keV") # Minimum energy
    emax = (11.0, "keV") # Maximum energy
    emission_field = "hard_emission" # The name of the field to use (normalization)
    alpha = 1.0 # The spectral index
    
    plaw_model = pyxsim.PowerLawSourceModel(e0, emin, emax, emission_field, alpha)
    
Another example where you have a spatially varying spectral index:

.. code-block:: python

    e0 = YTQuantity(2.0, "keV") # Reference energy
    emin = YTQuantity(0.2, "keV") # Minimum energy
    emax = YTQuantity(30.0, "keV") # Maximum energy
    emission_field = "inverse_compton_emission" # The name of the field to use (normalization)
    alpha = ("gas", "spectral_index") # The spectral index field
    
    plaw_model = pyxsim.PowerLawSourceModel(e0, emin, emax, emission_field, alpha)

.. _line-sources:

Line Emission Sources
---------------------

:class:`~pyxsim.source_models.LineSourceModel` assumes that the emission is 
occuring at a single energy, and that it may or may not be broadened by thermal
or other motions. In the former case, the emission is a delta function at a 
single rest-frame energy :math:`E_0`:

.. math::

    \varepsilon(E) = A\delta(E-E_0)

In the latter case, the emission is represented by a Gaussian with mean 
:math:`E_0` and standard deviation :math:`\sigma_E`:

.. math::

    \varepsilon(E) = \frac{A}{\sigma_E\sqrt{2\pi}}e^{-\frac{(E-E_0)^2}{2\sigma_E^2}}

When creating a :class:`~pyxsim.source_models.LineSourceModel`, it is 
initialized with the line rest-frame energy ``e0`` and an ``emission_field`` 
field specification that represents the normalization :math:`A` in the equations 
above, which must be in units of counts/s. Optionally, the line may be broadened 
by passing in a ``sigma`` parameter, which can be a field specification or 
``YTQuantity``, corresponding to either a spatially varying field or a single 
constant value. In either case, ``sigma`` may have units of energy or velocity;
if the latter, it will be converted to a broadening in energy units via 
:math:`\sigma_E = \sigma_v\frac{E_0}{c}`.

.. note:: 

    In most cases, you will want velocity broadening of lines to be handled by 
    the inputted velocity fields instead of by the ``sigma`` parameter. This 
    parameter is designed for thermal or other sources of "intrinsic" 
    broadening.

Examples
++++++++

An example of an unbroadened line:

.. code-block:: python

    e0 = YTQuantity(5.0, "keV") # Rest-frame line energy
    emission_field = ("gas", "line_emission") # Line emission field (normalization)
    line_model = pyxsim.LineSourceModel(e0, line_emission)

An example of a line with a constant broadening in km/s:

.. code-block:: python

    e0 = YTQuantity(6.0, "keV")
    emission_field = ("gas", "line_emission") # Line emission field (normalization)
    sigma = (500., "km/s")
    line_model = pyxsim.LineSourceModel(e0, line_emission, sigma=sigma)

An example of a line with a spatially varying broadening field:

.. code-block:: python

    e0 = YTQuantity(6.0, "keV")
    emission_field = ("gas", "line_emission") # Line emission field (normalization)
    sigma = "dark_matter_velocity_dispersion" # Has dimensions of velocity
    line_model = pyxsim.LineSourceModel(e0, line_emission, sigma=sigma)

Designing Your Own Source Model
-------------------------------

Though the three source models above cover a wide variety of possible use cases
for X-ray emission, you may find that you need to add a different source
altogether. It is possible to create your own source model to generate photon 
energies and positions. We will outline in brief the required steps to do so 
here. We'll use the already exising 
:class:`~pyxsim.source_models.PowerLawSourceModel` as an example.

To create a new source model, you'll need to make it a subclass of 
``SourceModel``. The first thing your source model needs is an ``__init__``
method to initialize a new instance of the model. This is where you pass in 
necessary parameters and initialize specific quantities such as the 
``spectral_norm`` and ``redshift`` to ``None``. These will be set to their 
appropriate values later, in the ``setup_model`` method. In this case, for 
a power-law spectrum, we need to define the maximum and minimum energies of the
spectrum (``emin`` and ``emax``), a reference energy (``e0``), an emissivity 
field that normalizes the spectrum (``emission_field``), and a spectral index 
field or single number ``alpha``:

.. code-block:: python

    def __init__(self, e0, emin, emax, emission_field, alpha, prng=None):
        self.e0 = parse_value(e0, "keV")
        self.emin = parse_value(emin, "keV")
        self.emax = parse_value(emax, "keV")
        self.emission_field = emission_field
        self.alpha = alpha
        self.prng = parse_prng(prng)
        self.spectral_norm = None
        self.redshift = None
        self.ftype = None

You need to also have an attribute for the yt field type stored in 
``self.ftype`` so that things such as position and velocity fields can be
determined. It's also always a good idea to have an optional keyword argument
``prng`` for a custom pseudo-random number generator. In this way, you can pass
in a random number generator (such as a :class:`~numpy.random.RandomState` 
instance) to get reproducible results. 

The next method you need to specify is the ``setup_model`` method:

.. code-block:: python

    def setup_model(self, data_source, redshift, spectral_norm):
        self.spectral_norm = spectral_norm
        self.redshift = redshift
        self.scale_factor = 1.0 / (1.0 + self.redshift)
        self.ftype = data_source.ds._get_field_info(self.emission_field).name[0]

It is called from :meth:`~pyxsim.photon_list.PhotonList.from_data_source` and is
used to set up the distance, redshift, and other aspects of the source being 
simulated. This does not happen in ``__init__`` because we may want to use the 
same source model for a number of different sources. You need to use one of the 
normalization fields (in this case the emission field) to determine the field
type.

The next method you need is ``__call__``. ``__call__`` is where the action 
really happens and the photon energies are generated. ``__call__`` takes a 
chunk of data from the data source, and for this chunk determines the emission
coming from each cell based on the normalization of the emission (in this case
given by the yt field ``"norm_field"``) and the spectrum of the source. We have
reproduced the method here with additional comments so that it is clearer
what is going on.

.. code-block:: python

    def __call__(self, chunk):

        # Determine the number of cells in this chunk
        num_cells = len(chunk[self.norm_field])

        # alpha can either be a single float number (the spectral index
        # is the same everywhere), or a spatially-dependent field.
        if isinstance(self.alpha, float):
            alpha = self.alpha*np.ones(num_cells)
        else:
            alpha = chunk[self.alpha].v

        # Here we are integrating the power-law spectrum over energy
        # between emin and emax. "norm_fac" represents the factor
        # you get when this is done. We need special logic here to
        # handle both the general case where alpha != 1 and where
        # alpha == 1. The "norm" that we compute at the end represents
        # the approximate number of photons in each cell.
        norm_fac = (self.emax.v**(1.-alpha)-self.emin.v**(1.-alpha))
        norm_fac[alpha == 1] = np.log(self.emax.v/self.emin.v)
        norm = norm_fac*chunk[self.emission_field].v*self.e0.v**alpha
        norm[alpha != 1] /= (1.-alpha[alpha != 1])
        norm *= self.spectral_norm*self.scale_factor

        # "norm" is now the approximate number of photons in each cell.
        # We will determine the number of photons from "norm" assuming
        # a Poisson distribution.
        number_of_photons = self.prng.poisson(lam=norm)

        # Generate an empty array for the energies
        energies = np.zeros(number_of_photons.sum())

        # Here we loop over the cells and determine the energies of the
        # photons in each cell by inverting the cumulative distribution
        # function corresponding to the power-law spectrum. Here again,
        # we have to do this differently depending on whether or not
        # alpha == 1.
        start_e = 0
        end_e = 0
        for i in range(num_cells):
            if number_of_photons[i] > 0:
                end_e = start_e+number_of_photons[i]
                u = self.prng.uniform(size=number_of_photons[i])
                if alpha[i] == 1:
                    e = self.emin.v*(self.emax.v/self.emin.v)**u
                else:
                    e = self.emin.v**(1.-alpha[i]) + u*norm_fac[i]
                    e **= 1./(1.-alpha[i])
                energies[start_e:end_e] = e * self.scale_factor
                start_e = end_e

        # Finally, __call__ must report the number of cells with photons, the 
        # number of photons in each cell which actually has photons, the actual 
        # indices of the cells themselves,
        # and the energies of the photons.
        active_cells = number_of_photons > 0
        ncells = active_cells.sum()

        return ncells, number_of_photons[active_cells], active_cells, energies[:end_e].copy()
