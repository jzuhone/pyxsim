Source Models for Generating Photons
====================================

pyXSIM comes with three pre-defined ``SourceModel`` types for 
generating photons. Though these should cover the vast majority of use cases, 
there is also the option to design your own source model. 

Thermal Sources
---------------

:class:`~pyxsim.source_models.ThermalSourceModel` assumes the emission of a hot 
thermal plasma can be described by a model that is only dependent on temperature 
and metallicity, and is proportional to the density squared:

.. math::

    \varepsilon(E) = n_en_H\Lambda(T, Z)


Power-Law Sources
-----------------

:class:`~pyxsim.source_models.PowerLawSourceModel` assumes that the emission can be 
described by a power law:

.. math::

    \varepsilon(E) = K\left(\frac{E}{E_0}\right)^{-\alpha}
    
between the energies ``emin`` and ``emax``, with a power-law spectral index ``alpha``,
normalized by an emission field specified by the user. 

Line Emission Sources
---------------------

:class:`~pyxsim.source_models.LineSourceModel` assumes that the emission is occuring at a 
single energy, and that it may or may not be broadened by thermal or other motions:

.. math::

    \varepsilon(E) = A\delta(E-E_0)

or:

.. math::

    \varepsilon(E) = \frac{A}{\sigma_E\sqrt{2\pi}}e^{-\frac{(E-E_0)^2}{2\sigma_E^2}}

The emissivity is normalized by an emission field specified by the user.

Designing Your Own Source Model
-------------------------------

We'll use the already exising :class:`~pyxsim.source_models.PowerLawSourceModel` as an example.

The first thing your source model needs is an ``__init__`` method to create an instance of the
model. This is where you pass in necessary parameters and initialize specific parameters such as
the ``spectral_norm`` and ``redshift`` to ``None``. These will be set to their appropriate values
later, in the ``setup_model`` method. In this case, we need to

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

The next method is ``__call__``:

.. code-block:: python

        def __call__(self, chunk):

            num_cells = len(chunk[self.norm_field])

            if isinstance(self.alpha, float):
                alpha = self.alpha*np.ones(num_cells)
            else:
                alpha = chunk[self.alpha].v

            norm_fac = (self.emax.v**(1.-alpha)-self.emin.v**(1.-alpha))
            norm_fac[alpha == 1] = np.log(self.emax.v/self.emin.v)
            norm = norm_fac*chunk[self.norm_field].v*self.e0.v**alpha
            norm[alpha != 1] /= (1.-alpha[alpha != 1])
            norm *= self.spectral_norm
            norm = np.modf(norm)

            u = self.prng.uniform(size=num_cells)
            number_of_photons = np.uint64(norm[1]) + np.uint64(norm[0] >= u)

            energies = np.zeros(number_of_photons.sum())

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
                    energies[start_e:end_e] = e / (1.+self.redshift)
                    start_e = end_e

            active_cells = number_of_photons > 0

            return number_of_photons[active_cells], active_cells, energies[:end_e].copy()

``__call__`` is where the action really happens and the photon energies are generated. In this case,
``__call__`` takes a chunk of data from the data source, and for this chunk determines the emission
coming from each cell based on the normalization of the emission (given by the yt field ``"norm_field"``)
and the spectral index of the source.

Finally, your source model needs a ``cleanup_model`` method to free memory, close file handles, and
reset the values of parameters that it used, in case you want to use the same source model instance
to generate photons for a different redshift, distance, etc. The ``cleanup_model`` method for
:class:`~pyxsim.source_models.PowerLawSourceModel` is very simple:

.. code-block:: python

        def cleanup_model(self):
            self.redshift = None
            self.spectral_norm = None
