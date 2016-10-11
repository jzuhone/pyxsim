.. _source-models:

Source Models for Generating Photons
====================================

pyXSIM comes with three pre-defined ``SourceModel`` types for generating a new
:class:`~pyxsim.photon_list.PhotonList`, for use with the 
:meth:`~pyxsim.photon_list.PhotonList.from_data_source` method. Though these 
should cover the vast majority of use cases, there is also the option to design
your own source model. 

.. _thermal-sources:

Thermal Sources
---------------

:class:`~pyxsim.source_models.ThermalSourceModel` assumes the emission of a hot 
thermal plasma can be described by a model that is only dependent on temperature 
and metallicity, and is proportional to the density squared:

.. math::

    \varepsilon(E) = n_en_H\Lambda(T, Z, E)

:class:`~pyxsim.source_models.ThermalSourceModel` requires the use of one of the
two thermal spectral models, described in the next sub-section. From this
spectral model, which depends on temperature and metallicity, a spectrum of photon
energies can be generated from each cell or particle. Since generating a new spectrum
for each cell would be prohibitively expensive, the cells are binned into a 1-D 
histogram of temperatures, and for each bin a spectrum is calculated. Provided the bins
are finely spaced enough, the accuracy of this method is sufficient for most purposes. 

There are a number of other optional parameters which can be set:

* ``temperature_field``: The yt field to use as the temperature. Must have dimensions
  of temperature. The default is ``("gas", "temperature")`` for grid-based datasets and
  ``("PartType0", "Temperature")`` or ``("io", "temperature")`` for particle-based datasets,
  depending on which is available.
* ``emission_measure_field``: The yt field to use as the emission measure. Must have
  dimensions of number density or per-volume. The default is ``("gas", "emission_measure")``
  for grid-based datasets. For particle-based datasets, a new field is constructed, using
  the default density and mass fields of the dataset, and the fields ``("PartType0", "ElectronAbundance")``
  ``("PartType0", "NeutralHydrogenAbundance")`` to construct the electron and hydrogen ion
  number densities if they are present in the dataset.
* ``kT_min``: The minimum temperature in units of keV in the set of temperature bins. Default is 0.0808.
* ``kT_max``: The maximum temperature in units of keV in the set of temperature bins. Default is 50.0.
* ``n_kT``: The number of temperature bins to use. Default is 10000.
* ``kT_scale``: The scaling of the temperature bins, either "linear" or "log". Default: "linear"
* ``Zmet``: The metallicity. Either a floating-point number for a constant metallicity, or the name of 
  a yt field for a spatially-varying metallicity. Default is 0.3.
* ``method``: The method used to generate the photon energies from the spectrum. Either ``"invert_cdf"``,
  which inverts the cumulative distribution function of the spectrum, or ``"accept_reject"``, which uses 
  the acceptance-rejection method on the spectrum. The first method should be sufficient for most cases. 
* ``prng``: A pseudo-random number generator. Typically will only be specified
  if you have a reason to generate the same set of random numbers, such as for a 
  test or a comparison. Default is the :mod:`numpy.random` module, but a 
  :class:`~numpy.random.RandomState` object can also be used. 

Thermal Spectra
+++++++++++++++

For thermal plasmas, pyXSIM provides two classes for modeling emission that is dependent on temperature
and metallicity, e.g. :math:`\Lambda(T,Z)`. These can serve as input to the 
:class:`~pyxsim.source_models.ThermalSourceModel` class. Both of these classes are 
:class:`~pyxsim.spectral_models.ThermalSpectralModel` objects.

:class:`~pyxsim.spectral_models.TableApecModel` generates a thermal emission spectrum
from the APEC plasma emission tables available from `AtomDB <http://www.atomdb.org>`_:

.. code-block:: python

    emin = 0.01 # The minimum energy of the spectrum in keV
    emax = 20.0 # The maximum energy of the spectrum in keV
    nchan = 10000 # The number of spectral channels
    apec_vers = "3.0.3" # The version identifier string for the APEC files. Default: "2.0.2"
    apec_root = "/Users/jzuhone/atomdb" # The directory where the APEC model files are stored.
                                        # Optional, the native pyxsim files will be used if
                                        # a location is not provided.
    spec_model = pyxsim.TableApecModel(emin, emax, nchan, apec_root=apec_root,
                                       apec_vers=apec_vers, thermal_broad=False)

:class:`~pyxsim.spectral_models.XSpecThermalModel` generates a thermal emission spectrum
from a model known to XSPEC, using PyXspec as a backend:

.. code-block:: python

    model = "apec" # or "mekal", "bapec", etc.
    emin = 0.01 # The minimum energy of the spectrum in keV
    emax = 20.0 # The maximum energy of the spectrum in keV
    nchan = 10000 # The number of spectral channels
    spec_model = pyxsim.XSpecThermalModel(model, emin, emax, nchan, thermal_broad=True)

They keyword argument ``thermal_broad`` should be set to ``True`` or ``False`` depending on
whether or not you want the spectral lines thermally broadened. 

You will need to set up one of these two models in your script and pass it as the first argument to
:class:`~pyxsim.source_models.ThermalSourceModel`.

Though a :class:`~pyxsim.spectral_models.ThermalSpectralModel` is mainly used internally by the 
:class:`~pyxsim.source_models.ThermalSourceModel` to construct spectra, there is also a method
:meth:`~pyxsim.source_models.ThermalSourceModel.return_spectrum` which can be used to return a 
spectrum (in an array) for a given temperature, metallicity, redshift, and flux normalization:

.. code-block:: python

    kT = 6.0 # in keV
    metallicity = 0.3 # in solar units
    z = 0.05 # redshift
    norm = 1.0e-3 # total spectrum flux in units of photons/s/cm**2
    
    spec = spec_model.return_spectrum(kT, metallicity, z, norm)

The units of the returned spectrum are in :math:`{\rm photons~s^{-1}~cm^{-1}}`.

Examples
++++++++

Here, we will show several examples of constructing :class:`~pyxsim.source_models.ThermalSourceModel`
objects. First, construct a thermal spectral model:

.. code-block:: python

    spec_model = pyxsim.TableApecModel(0.1, 20.0, 10000, thermal_broad=True)

An example where we use the default parameters, except we set a constant metallicity:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel(spec_model, Zmet=0.5)

An example where we use a metallicity field and change the temperature field:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel(spec_model, Zmet=("gas", "metallicity"),
                                              temperature_field="hot_gas_temp")

An example where we change the limits and number of the temperature bins:

.. code-block:: python

    thermal_model = pyxsim.ThermalSourceModel(spec_model, kT_min=0.1, kT_max=100.,
                                              n_kT=50000)

An example where we specify a random number generator:

.. code-block:: python

    from numpy.random import RandomState
    prng = RandomState(25)
    thermal_model = pyxsim.ThermalSourceModel(spec_model, prng=prng)

.. _power-law-sources:

Power-Law Sources
-----------------

:class:`~pyxsim.source_models.PowerLawSourceModel` assumes that the emission can be 
described by a pure power law:

.. math::

    \varepsilon(E) = K\left(\frac{E}{E_0}\right)^{-\alpha}, E_{\rm min} \leq E \leq E_{\rm max}
    
between the energies ``emin`` and ``emax``, with a power-law spectral index ``alpha``.
The power law normalization :math:`K` is represented by an ``emission_field`` specified 
by the user, which must have units of counts/s/keV in the source rest frame. ``alpha``
may be a single floating-point number (implying the spectral index is the same everywhere), 
or a field specification corresponding to a spatially varying spectral index. A reference
energy ``e0`` (see above equation) must also be specified.

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

:class:`~pyxsim.source_models.LineSourceModel` assumes that the emission is occuring at a 
single energy, and that it may or may not be broadened by thermal or other motions. In the 
former case, the emission is a delta function at a single rest-frame energy :math:`E_0`:

.. math::

    \varepsilon(E) = A\delta(E-E_0)

In the latter case, the emission is represented by a Gaussian with mean :math:`E_0` and
standard deviation :math:`\sigma_E`:

.. math::

    \varepsilon(E) = \frac{A}{\sigma_E\sqrt{2\pi}}e^{-\frac{(E-E_0)^2}{2\sigma_E^2}}

When creating a :class:`~pyxsim.source_models.LineSourceModel`, it is initialized with
the line rest-frame energy ``e0`` and an ``emission_field`` field specification that 
represents the normalization :math:`A` in the equations above, which must be in units of
counts/s. Optionally, the line may be broadened by passing in a ``sigma`` parameter, which
can be a field specification or ``YTQuantity``, corresponding to either a spatially
varying field or a single constant value. In either case, ``sigma`` may have units of energy or 
velocity; if the latter, it will be converted to a broadening in energy units via
:math:`\sigma_E = \sigma_v\frac{E_0}{c}`.

.. note:: 

    In most cases, you will want velocity broadening of lines to be handled by the 
    inputted velocity fields instead of by the ``sigma`` parameter. This parameter
    is designed for thermal or other sources of "intrinsic" broadening.

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

Though the three source models above cover a wide variety of possible use cases for X-ray emission,
you may find that you need to add a different source altogether. It is possible to create your own
source model to generate photon energies and positions. We will outline in brief the required steps
to do so here. We'll use the already exising :class:`~pyxsim.source_models.PowerLawSourceModel` as
an example.

To create a new source model, you'll need to make it a subclass of ``SourceModel``. The first thing
your source model needs is an ``__init__`` method to initialize a new instance of the model. This is
where you pass in necessary parameters and initialize specific quantities such as the ``spectral_norm``
and ``redshift`` to ``None``. These will be set to their appropriate values later, in the ``setup_model``
method. In this case, for a power-law spectrum, we need to define the maximum and minimum energies of the
spectrum (``emin`` and ``emax``), a reference energy (``e0``), an emissivity field that normalizes the
spectrum (``norm_field``), and a spectral index field or single number ``alpha``:

.. code-block:: python

    class PowerLawSourceModel(SourceModel):
        def __init__(self, e0, emin, emax, norm_field, alpha, prng=None):
            self.e0 = parse_value(e0, "keV")
            self.emin = parse_value(emin, "keV")
            self.emax = parse_value(emax, "keV")
            self.norm_field = norm_field
            self.alpha = alpha
            if prng is None:
                self.prng = np.random
            else:
                self.prng = prng
            self.spectral_norm = None
            self.redshift = None

It's also always a good idea to have an optional keyword argument ``prng`` for a custom pseudo-random
number generator. In this way, you can pass in a random number generator (such as a :class:`~numpy.random.RandomState`
instance) to get reproducible results. The default should be the :mod:`~numpy.random` module.

The next method you need to specify is the ``setup_model`` method:

.. code-block:: python

    def setup_model(self, data_source, redshift, spectral_norm):
        self.spectral_norm = spectral_norm
        self.redshift = redshift

``setup_model`` should always have this exact method signature. It is called from :meth:`~pyxsim.photon_list.PhotonList.from_data_source`
and is used to set up the distance, redshift, and other aspects of the source being simulated. This does not happen in
``__init__`` because we may want to use the same source model for a number of different sources.

The next method you need is ``__call__``. ``__call__`` is where the action really happens and the photon energies
are generated. ``__call__`` takes a chunk of data from the data source, and for this chunk determines the emission
coming from each cell based on the normalization of the emission (in this case given by the yt field ``"norm_field"``)
and the spectrum of the source. We have reproduced the method here with additional comments so that it is clearer
what is going on.

.. code-block:: python

    def __call__(self, chunk):

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
        norm = norm_fac*chunk[self.norm_field].v*self.e0.v**alpha
        norm[alpha != 1] /= (1.-alpha[alpha != 1])
        norm *= self.spectral_norm

        # "norm" is now the approximate number of photons in each cell.
        # what we have to do next is determine the actual number of
        # photons in each cell. What we do here is split "norm" into
        # its integer and fractional parts, and use the latter as the
        # probability that an extra photon will be observed from this
        # cell in addition to those from the integer part.
        norm = np.modf(norm)
        u = self.prng.uniform(size=num_cells)
        number_of_photons = np.uint64(norm[1]) + np.uint64(norm[0] >= u)

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
                # Scale by the redshift
                energies[start_e:end_e] = e / (1.+self.redshift)
                start_e = end_e

        # Finally, __call__ must report the number of photons in each cell
        # which actually has photons, the actual indices of the cells themselves,
        # and the energies of the photons.
        active_cells = number_of_photons > 0

        return number_of_photons[active_cells], active_cells, energies[:end_e].copy()

Finally, your source model needs a ``cleanup_model`` method to free memory, close file handles, and
reset the values of parameters that it used, in case you want to use the same source model instance
to generate photons for a different redshift, distance, etc. The ``cleanup_model`` method for
:class:`~pyxsim.source_models.PowerLawSourceModel` is very simple:

.. code-block:: python

    def cleanup_model(self):
        self.redshift = None
        self.spectral_norm = None
