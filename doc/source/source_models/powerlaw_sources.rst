.. _powerlaw-sources:

Power-Law Sources
-----------------

:class:`~pyxsim.source_models.PowerLawSourceModel` assumes that the emission can
be described by a pure power law:

.. math::

    \varepsilon(E) = K\left(\frac{E}{E_0}\right)^{-\alpha}, E_{\rm min} \leq E \leq E_{\rm max}

between the energies ``emin`` and ``emax``, with a power-law spectral index
``alpha``. The power law normalization K is represented by a ``luminosity_field``
specified by the user, which represents the per-cell or per-particle luminosity
in the source rest frame within the ``emin``-``emax`` band, and must have units
of power, such as W, erg/s, or keV/s. ``alpha`` may be a single floating-point
number (implying the spectral index is the same everywhere), or a field
specification corresponding to a spatially varying spectral index. A reference
energy ``e0`` (see above equation) must also be specified.

Examples
++++++++

An example where the spectral index is the same everywhere:

.. code-block:: python

    e0 = (1.0, "keV") # Reference energy
    emin = (0.01, "keV") # Minimum energy
    emax = (11.0, "keV") # Maximum energy
    luminosity_field = "hard_emission" # The name of the field to use (normalization)
    alpha = 1.0 # The spectral index

    plaw_model = pyxsim.PowerLawSourceModel(e0, emin, emax, luminosity_field, alpha)

Another example where you have a spatially varying spectral index:

.. code-block:: python

    e0 = YTQuantity(2.0, "keV") # Reference energy
    emin = YTQuantity(0.2, "keV") # Minimum energy
    emax = YTQuantity(30.0, "keV") # Maximum energy
    luminosity_field = "inverse_compton_luminosity" # The name of the field to use (normalization)
    alpha = ("gas", "spectral_index") # The spectral index field

    plaw_model = pyxsim.PowerLawSourceModel(e0, emin, emax, luminosity_field, alpha)
