from pyxsim.spectral_models import CXSpectralModel

from .base import ThermalSourceModel


class CXSourceModel(ThermalSourceModel):
    """
    This class generates a source model for emission from charge exchange,
    using the AtomDB Charge Exchange Model, v2.0 (ACX2). It models the
    emission obtained from charges exchanged between donor neutral hydrogen/helium
    atoms colliding with recombining ions. The hydrogen number density for the
    recombining plasma, and the fraction of neutral helium must be supplied, as
    detailed below. The "collision parameter" must also be specified, which is the
    relative velocity between the ions and the neutrals. This can be either
    a single value or a spatially varying field. The emission spectrum is a
    function of this parameter and is interpolated from a precomputed table
    for each ion, the velocity bins for which can also be specified below.
    Other ACX2 parameters can also be set. For the meaning of these parameters,
    consult the ACX2 docs at https://acx2.readthedocs.io/.

    To use this model, you must have the AtomDB Charge Exchange Model
    package installed from https://github.com/AtomDB/ACX2, as well as the
    pyatomdb package.

    Parameters
    ----------
    emin : float
        The minimum energy for the spectrum in keV.
    emax : float
        The maximum energy for the spectrum in keV.
    nbins : integer
        The number of channels in the spectrum.
    collnpar : float, (value, unit) tuple, unyt_quantity, or Quantity
        The collision parameter for the CX process, in units of velocity.
        If a float, the units are assumed to be km/s. If set to a field name,
        this will be a spatially varying collision parameter field.
    h_r_number_density : tuple of strings
        The field representing the number density of hydrogen in the
        recombining gas.
    he_d_fraction : float, string, or tuple of strings, optional
        The neutral helium mass fraction in the donor gas. If a float, this
        assumes a constant mass fraction throughout. If a tuple of strings, it
        is taken to be the name of the given field representing the neutral
        helium fraction. This assumes that the total helium and hydrogen mass
        fractions in the donor gas sum to 1.
    Zmet : float, string, or tuple of strings
        The metallicity. If a float, assumes a constant metallicity throughout
        in solar units. If a string or tuple of strings, is taken to be the
        name of the metallicity field.
    var_elem : dict, optional
        Elements that should be allowed to vary freely from the single abundance
        parameter. Each dictionary value, specified by the abundance symbol,
        corresponds to the abundance of that symbol. If a float, it is understood
        to be constant and in solar units. If a string or tuple of strings, it is
        assumed to be a spatially varying field. Default: None
    acx_model : integer, optional
        ACX model to fall back on, from 1 to 8. Default: 8.
    recomb_type : integer, optional
        Single recombination (1) or all the way to neutral (2) Default: 1.
    vmin : float, optional
        The minimum value of the velocity table in km/s. Default: 10.0
    vmax : float, optional
        The maximum value of the velocity table in km/s. Default: 10000.0
    nbins_v : integer, optional
        The number of bins in the velocity table. Default: 100
    binscale : string, optional
        The scale of the energy binning: "linear" or "log".
        Default: "linear"
    temperature_field : string or (ftype, fname) tuple, optional
        The yt temperature field to use for the thermal modeling. Must have
        units of Kelvin. Default: ("gas","temperature")
    h_fraction : float, string, or tuple of strings, optional
        The total hydrogen mass fraction. If a float, assumes a constant mass
        fraction of hydrogen throughout. If a string or tuple of strings,
        is taken to be the name of the hydrogen fraction field. Default is
        whatever value is appropriate for the chosen abundance table.
    kT_min : float, optional
        The default minimum temperature in keV to compute emission for.
        Default: 0.025
    kT_max : float, optional
        The default maximum temperature in keV to compute emission for.
        Default: 64.0
    max_density : float, (value, unit) tuple, unyt_quantity, or Quantity
        The maximum density of the cells or particles to use when generating
        photons. If a float, the units are assumed to be g/cm**3.
        Default: None, meaning no maximum density.
    min_entropy : float, (value, unit) tuple, unyt_quantity, or Quantity
        The minimum entropy of the cells or particles to use when generating
        photons. If a float, the units are assumed to be keV*cm**2.
        Default: None, meaning no minimum entropy.
    method : string, optional
        The method used to generate the photon energies from the spectrum:
        "invert_cdf": Invert the cumulative distribution function of the spectrum.
        "accept_reject": Acceptance-rejection method using the spectrum.
        The first method should be sufficient for most cases.
    abund_table : string or array_like, optional
        The abundance table to be used for solar abundances.
        Either a string corresponding to a built-in table or an array
        of 30 floats corresponding to the abundances of each element
        relative to the abundance of H. Default is "angr".
        Built-in options are:
        "angr" : from Anders E. & Grevesse N. (1989, Geochimica et
        Cosmochimica Acta 53, 197)
        "aspl" : from Asplund M., Grevesse N., Sauval A.J. & Scott
        P. (2009, ARAA, 47, 481)
        "wilm" : from Wilms, Allen & McCray (2000, ApJ 542, 914
        except for elements not listed which are given zero abundance)
        "lodd" : from Lodders, K (2003, ApJ 591, 1220)
        "feld" : from Feldman U. (Physica Scripta, 46, 202)
        "cl17.03" : the abundance table used in Cloudy v17.03.
    prng : integer or numpy.random.RandomState object
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is to use the :mod:`numpy.random` module.

    Examples
    --------
    >>> source_model = CXSourceModel(
    ...     0.1, 3.0, 10000, ("gas", "coll_vel"), ("gas", "h_p1_number_density"),
    ...     ("gas", "he_p0_fraction"), Zmet)
    """

    _cx = True

    def __init__(
        self,
        emin,
        emax,
        nbins,
        collnpar,
        h_r_number_density,
        h_d_number_density,
        he_d_number_density,
        Zmet,
        var_elem=None,
        acx_model=8,
        recomb_type=1,
        vmin=10.0,
        vmax=10000.0,
        nbins_v=100,
        binscale="linear",
        temperature_field=("gas", "temperature"),
        h_fraction=None,
        kT_min=0.01,
        kT_max=64.0,
        nbins_kT=100,
        max_density=None,
        min_entropy=None,
        method="invert_cdf",
        abund_table="angr",
        prng=None,
    ):
        spectral_model = CXSpectralModel(
            emin,
            emax,
            nbins,
            vmin,
            vmax,
            nbins_v,
            kT_min,
            kT_max,
            nbins_kT,
            collntype=2,
            acx_model=acx_model,
            recomb_type=recomb_type,
            binscale=binscale,
            abund_table=abund_table,
            var_elem=var_elem,
        )
        super().__init__(
            spectral_model,
            emin,
            emax,
            nbins,
            Zmet,
            binscale=binscale,
            kT_min=kT_min,
            kT_max=kT_max,
            var_elem=var_elem,
            max_density=max_density,
            min_entropy=min_entropy,
            method=method,
            abund_table=abund_table,
            prng=prng,
            h_fraction=h_fraction,
            temperature_field=temperature_field,
        )
        self.collnpar = collnpar
        self.h_r_number_density = h_r_number_density
        self.h_d_number_density = h_d_number_density
        self.he_d_number_density = he_d_number_density
        self.collntype = 2
        self.acx_model = acx_model
        self.recomb_type = recomb_type
        self.var_elem_keys = self.spectral_model.var_elem_names
        self.var_ion_keys = self.spectral_model.var_ion_names

    def _prep_repr(self):
        class_name, strs = super()._prep_repr()
        strs["collnpar"] = self.collnpar
        strs["collntype"] = self.collntype
        strs["acx_model"] = self.acx_model
        strs["recomb_type"] = self.recomb_type
        strs.pop("model")
        return class_name, strs


class CXNEISourceModel(CXSourceModel):
    """
    This class is almost identical to the CXSourceModel class, except that
    one must specify the abundances of the various elemental ion states by
    hand, similar to the NEISourceModel.

    This class generates a source model for emission from charge exchange,
    using the AtomDB Charge Exchange Model, v2.0 (ACX2). It models the
    emission obtained from charges exchanged between donor neutral hydrogen/helium
    atoms colliding with recombining ions. The hydrogen number density for the
    recombining plasma, and the fraction of neutral helium must be supplied, as
    detailed below. The "collision parameter" must also be specified, which is the
    relative velocity between the ions and the neutrals. This can be either
    a single value or a spatially varying field. The emission spectrum is a
    function of this parameter and is interpolated from a precomputed table
    for each ion, the velocity bins for which can also be specified below.
    Other ACX2 parameters can also be set. For the meaning of these parameters,
    consult the ACX2 docs at https://acx2.readthedocs.io/.

    To use this model, you must have the AtomDB Charge Exchange Model
    package installed from https://github.com/AtomDB/ACX2, as well as the
    pyatomdb package.

    Parameters
    ----------
    emin : float
        The minimum energy for the spectrum in keV.
    emax : float
        The maximum energy for the spectrum in keV.
    nbins : integer
        The number of channels in the spectrum.
    collnpar : float, (value, unit) tuple, unyt_quantity, or Quantity
        The collision parameter for the CX process, in units of velocity.
        If a float, the units are assumed to be km/s. If set to a field name,
        this will be a spatially varying collision parameter field.
    h_r_number_density : tuple of strings
        The field representing the number density of hydrogen in the
        recombining gas.
    he_d_fraction : float, string, or tuple of strings, optional
        The neutral helium mass fraction in the donor gas. If a float, this
        assumes a constant mass fraction throughout. If a tuple of strings, it
        is taken to be the name of the given field representing the neutral
        helium fraction. This assumes that the total helium and hydrogen mass
        fractions in the donor gas sum to 1.
    var_elem : dict
        Abundances of ions. Each dictionary value, specified by the ionic
        symbol, corresponds to the abundance of that symbol. If a float, it is
        understood to be constant and in solar units. If a string or tuple of
        strings, it is assumed to be a spatially varying field.
    acx_model : integer, optional
        ACX model to fall back on, from 1 to 8. Default: 8.
    recomb_type : integer, optional
        Single recombination (1) or all the way to neutral (2) Default: 1.
    vmin : float, optional
        The minimum value of the velocity table in km/s. Default: 10.0
    vmax : float, optional
        The maximum value of the velocity table in km/s. Default: 10000.0
    nbins_v : integer, optional
        The number of bins in the velocity table. Default: 100
    binscale : string, optional
        The scale of the energy binning: "linear" or "log".
        Default: "linear"
    temperature_field : string or (ftype, fname) tuple, optional
        The yt temperature field to use for the thermal modeling. Must have
        units of Kelvin. Default: ("gas","temperature")
    h_fraction : float, string, or tuple of strings, optional
        The total hydrogen mass fraction. If a float, assumes a constant mass
        fraction of hydrogen throughout. If a string or tuple of strings,
        is taken to be the name of the hydrogen fraction field. Default is
        whatever value is appropriate for the chosen abundance table.
    kT_min : float, optional
        The default minimum temperature in keV to compute emission for.
        Default: 0.025
    kT_max : float, optional
        The default maximum temperature in keV to compute emission for.
        Default: 64.0
    max_density : float, (value, unit) tuple, unyt_quantity, or Quantity
        The maximum density of the cells or particles to use when generating
        photons. If a float, the units are assumed to be g/cm**3.
        Default: None, meaning no maximum density.
    min_entropy : float, (value, unit) tuple, unyt_quantity, or Quantity
        The minimum entropy of the cells or particles to use when generating
        photons. If a float, the units are assumed to be keV*cm**2.
        Default: None, meaning no minimum entropy.
    method : string, optional
        The method used to generate the photon energies from the spectrum:
        "invert_cdf": Invert the cumulative distribution function of the spectrum.
        "accept_reject": Acceptance-rejection method using the spectrum.
        The first method should be sufficient for most cases.
    abund_table : string or array_like, optional
        The abundance table to be used for solar abundances.
        Either a string corresponding to a built-in table or an array
        of 30 floats corresponding to the abundances of each element
        relative to the abundance of H. Default is "angr".
        Built-in options are:
        "angr" : from Anders E. & Grevesse N. (1989, Geochimica et
        Cosmochimica Acta 53, 197)
        "aspl" : from Asplund M., Grevesse N., Sauval A.J. & Scott
        P. (2009, ARAA, 47, 481)
        "wilm" : from Wilms, Allen & McCray (2000, ApJ 542, 914
        except for elements not listed which are given zero abundance)
        "lodd" : from Lodders, K (2003, ApJ 591, 1220)
        "feld" : from Feldman U. (Physica Scripta, 46, 202)
        "cl17.03" : the abundance table used in Cloudy v17.03.
    prng : integer or numpy.random.RandomState object
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is to use the :mod:`numpy.random` module.

    Examples
    --------
    >>> var_elem = {
    ...     "Si^12": ("flash", "si12"),
    ...     "S^14":  ("flash", "s14 "),
    ...     "Ar^16": ("flash", "ar16"),
    ...     "Ca^18": ("flash", "ca18"),
    ... }
    >>> source_model = CXNEISourceModel(
    ...     0.1, 3.0, 10000, ("gas", "coll_vel"), ("gas", "h_p1_number_density"),
    ...     ("gas", "he_p0_fraction"), var_elem)
    """

    _nei = True

    def __init__(
        self,
        emin,
        emax,
        nbins,
        collnpar,
        h_r_number_density,
        h_d_number_density,
        he_d_number_density,
        var_elem,
        acx_model=8,
        recomb_type=1,
        vmin=10.0,
        vmax=10000.0,
        nbins_v=100,
        binscale="linear",
        temperature_field=("gas", "temperature"),
        h_fraction=None,
        kT_min=0.01,
        kT_max=64.0,
        max_density=None,
        min_entropy=None,
        method="invert_cdf",
        abund_table="angr",
        prng=None,
    ):
        super().__init__(
            emin,
            emax,
            nbins,
            collnpar,
            h_r_number_density,
            h_d_number_density,
            he_d_number_density,
            0.0,
            var_elem=var_elem,
            acx_model=acx_model,
            recomb_type=recomb_type,
            vmin=vmin,
            vmax=vmax,
            nbins_v=nbins_v,
            binscale=binscale,
            temperature_field=temperature_field,
            h_fraction=h_fraction,
            kT_min=kT_min,
            kT_max=kT_max,
            max_density=max_density,
            min_entropy=min_entropy,
            method=method,
            abund_table=abund_table,
            prng=prng,
        )

    def _prep_repr(self):
        class_name, strs = super()._prep_repr()
        strs.pop("Zmet")
        return class_name, strs
