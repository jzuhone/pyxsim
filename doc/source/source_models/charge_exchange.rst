.. _charge-exchange:

Charge Exchange Sources
-----------------------

Charge exchange is the process where ions collide with neutral hydrogen and/or
helium atoms, transferring an electron to the neutral atom, often in an excited
state, which then decays to a lower state, producing line emission.

pyXSIM provides a model for charge exchange emission using the AtomDB Charge
Exchange Model, v2.0 (known as "ACX2" for short), which models this emission
obtained from collisions between neutral hydrogen/helium and ions. As with the
other pyXSIM source models, one needs fields to define the various elements and
physical quantities involved in the interaction.

    ions. The neutrals and the ion fields must be supplied, as detailed
    below. The "collision parameter" must also be specified, which is the
    relative velocity between the ions and the neutrals. This can be either
    a single value or a spatially varying field. The emission spectrum is a
    function of this parameter and is interpolated from a precomputed table
    for each ion, the velocity bins for which can also be specified below.
    Other ACX2 parameters can also be set. For the meaning of these parameters,
    consult the ACX2 docs at https://acx2.readthedocs.io/.

To use this model, you must have the
`AtomDB Charge Exchange Model <https://acx2.readthedocs.io/>`_ package
installed, as well as the `pyatomdb <https://pyatomdb.readthedocs.io>`_ package.
