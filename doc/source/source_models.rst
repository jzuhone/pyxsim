Source Models for Generating Photons
====================================

pyXSIM comes with three pre-defined ``SourceModel`` types for 
generating photons. Though these should cover the vast majority of use cases, 
there is also the option to design your own source model. 

Thermal Sources
---------------

``ThermalSourceModel`` assumes the emission of a hot thermal plasma can be described by a
model that is only dependent on temperature and metallicity, and is proportional
to the density squared:

.. math::

    \varepsilon(E) = n_en_H\Lambda(T, Z)


Power-Law Sources
-----------------

``PowerLawSourceModel`` assumes that the emission can be described by a power law:

.. math::

    \varepsilon(E) = K\left(\frac{E}{E_0}\right)^{-\alpha}
    
between the energies ``emin`` and ``emax``, with a power-law spectral index ``alpha``,
normalized by an emission field specified by the user. 

Line Emission Sources
---------------------

``LineSourceModel`` assumes that the emission is occuring at a single energy, and that
it may or may not be broadened by thermal or other motions:

.. math::

    \varepsilon(E) = A\delta(E-E_0)

or:

.. math::

    \varepsilon(E) = \frac{A}{\sigma_E\sqrt{2\pi}}e^{-\frac{(E-E_0)^2}{2\sigma_E^2}}

The emissivity is normalized by an emission field specified by the user. 