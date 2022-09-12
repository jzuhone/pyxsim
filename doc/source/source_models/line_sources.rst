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
