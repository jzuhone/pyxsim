from pyxsim.spectral_models import CloudyCIESpectralModel, MekalSpectralModel, TableCIEModel
from pyxsim.utils import mylog

from .base import ThermalSourceModel


class CIESourceModel(ThermalSourceModel):
    """
    Initialize a source model from a CIE spectrum, using either
    the APEC, SPEX, MeKaL, or Cloudy models.

    Parameters
    ----------
    model : string
        Which spectral emission model to use. Accepts either "apec", "spex",
        "mekal", or "cloudy".
    emin : float
        The minimum energy for the spectrum in keV.
    emax : float
        The maximum energy for the spectrum in keV.
    nbins : integer
        The number of channels in the spectrum.
    Zmet : float, string, or tuple of strings
        The metallicity. If a float, assumes a constant metallicity throughout
        in solar units. If a string or tuple of strings, is taken to be the
        name of the metallicity field.
    binscale : string, optional
        The scale of the energy binning: "linear" or "log".
        Default: "linear"
    temperature_field : string or (ftype, fname) tuple, optional
        The yt temperature field to use for the thermal modeling. Must have
        units of Kelvin. Default: ("gas", "temperature")
    emission_measure_field : str or (ftype, fname) tuple, optional
        The yt emission measure field to use for the thermal modeling. Must
        have units of cm^-3. Default: ("gas", "emission_measure")
    h_fraction : float, str, or tuple of strings, optional
        The hydrogen mass fraction. If a float, assumes a constant mass
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
    var_elem : dict, optional
        Elements that should be allowed to vary freely from the single abundance
        parameter. Each dictionary value, specified by the abundance symbol,
        corresponds to the abundance of that symbol. If a float, it is understood
        to be constant and in solar units. If a string or tuple of strings, it is
        assumed to be a spatially varying field. Not yet available for "cloudy".
        Default: None
    method : string, optional
        The method used to generate the photon energies from the spectrum:
        "invert_cdf": Invert the cumulative distribution function of the spectrum.
        "accept_reject": Acceptance-rejection method using the spectrum.
        The first method should be sufficient for most cases.
    thermal_broad : boolean, optional
        Whether the spectral lines should be thermally
        broadened. Only available for "apec" or "spex". Default: True
    model_root : string, optional
        The directory root where the model files are stored. If not provided,
        a default location known to pyXSIM is used.
    model_vers : string, optional
        The version identifier string for the model files, e.g.
        "2.0.2", if supported by the model. Currently only supported by
        "apec", "spex", and "cloudy". Default depends on the model being used.
        If "cloudy", the options are:
        "4_lo": Tables computed from Cloudy using a continuum resolution
        of 0.1 with a range of 0.05 to 10 keV.
        "4_hi": Tables computed from Cloudy using enhanced continuum
        resolution of 0.025 with a range of 0.05 to 10 keV. Excellent
        energy resolution, but may be expensive to evaluate. Default for
        "cloudy" is "4_lo".
    nolines : boolean, optional
        Turn off lines entirely for generating emission. Only available
        for "apec" or "spex". Default: False
    abund_table : string or array_like, optional
        The abundance table to be used for solar abundances.
        Either a string corresponding to a built-in table or an array
        of 30 floats corresponding to the abundances of each element
        relative to the abundance of H. Not yet available for the "cloudy",
        CIE model, which always uses "feld". Otherwise, default is "angr".
        Built-in options are:
        "angr" : from Anders E. & Grevesse N. (1989, Geochimica et
        Cosmochimica Acta 53, 197)
        "aspl" : from Asplund M., Grevesse N., Sauval A.J. & Scott
        P. (2009, ARAA, 47, 481)
        "wilm" : from Wilms, Allen & McCray (2000, ApJ 542, 914
        except for elements not listed which are given zero abundance)
        "lodd" : from Lodders, K (2003, ApJ 591, 1220)
        "feld" : from Feldman U. (1992, Physica Scripta, 46, 202)
        "cl17.03" : the abundance table used in Cloudy v17.03.
    trace_abund : float, optional
        The abundance to give to trace elements (Li, Be, B, F, Na, P,
        Cl, K, Sc, Ti, V, Cr, Mn, Co, Cu, Zn), relative to solar. Any
        trace element that has an abundance already set using var_elem
        will not be considered here. By default, trace element abundances
        are set at 1 solar, similar to the behavior of XSPEC. This option
        is only available for the "apec" and "spex" models.
    prng : integer or :class:`~numpy.random.RandomState` object
        A pseudo-random number generator. Typically, will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is to use the :mod:`numpy.random` module.

    Examples
    --------
    >>> source_model = CIESourceModel("apec", 0.1, 10.0, 10000,
    ...                               ("gas", "metallicity"))
    """

    _nei = False
    _density_dependence = False

    def __init__(
        self,
        model,
        emin,
        emax,
        nbins,
        Zmet,
        binscale="linear",
        temperature_field=("gas", "temperature"),
        emission_measure_field=("gas", "emission_measure"),
        h_fraction=None,
        kT_min=0.025,
        kT_max=64.0,
        max_density=None,
        min_entropy=None,
        var_elem=None,
        method="invert_cdf",
        thermal_broad=True,
        model_root=None,
        model_vers=None,
        nolines=False,
        abund_table="angr",
        trace_abund=1.0,
        prng=None,
    ):
        var_elem_keys = list(var_elem.keys()) if var_elem else None
        if model in ["apec", "spex"]:
            spectral_model = TableCIEModel(
                model,
                emin,
                emax,
                nbins,
                kT_min,
                kT_max,
                binscale=binscale,
                var_elem=var_elem_keys,
                thermal_broad=thermal_broad,
                model_root=model_root,
                model_vers=model_vers,
                nolines=nolines,
                nei=self._nei,
                abund_table=abund_table,
                trace_abund=trace_abund,
            )
            self.trace_abund = trace_abund
        elif model == "mekal":
            spectral_model = MekalSpectralModel(emin, emax, nbins, binscale=binscale, var_elem=var_elem_keys)
        elif model == "cloudy":
            if abund_table != "feld":
                mylog.warning(
                    "For the 'cloudy' model, the only available abundance table is 'feld', so using that one."
                )
                abund_table = "feld"
            spectral_model = CloudyCIESpectralModel(
                emin,
                emax,
                nbins,
                binscale=binscale,
                var_elem=var_elem_keys,
                model_vers=model_vers,
            )
        self.model = model
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
            temperature_field=temperature_field,
            emission_measure_field=emission_measure_field,
            h_fraction=h_fraction,
        )
        self.var_elem_keys = self.spectral_model.var_elem_names
        self.var_ion_keys = self.spectral_model.var_ion_names
        self.nolines = nolines
        self.thermal_broad = thermal_broad

    def _prep_repr(self):
        class_name, super_strs = super()._prep_repr()
        strs = {"model": self.model}
        strs.update(super_strs)
        strs["nolines"] = self.nolines
        strs["thermal_broad"] = self.thermal_broad
        strs["trace_abund"] = self.trace_abund
        return class_name, strs


class NEISourceModel(CIESourceModel):
    """
    Initialize a source model from a thermal spectrum, using the
    APEC NEI tables from https://www.atomdb.org. Note that for this
    class a set of specific element fields must be supplied. This should
    really only be used with simulation data which include the appropriate
    NEI calculations.

    Parameters
    ----------
    emin : float
        The minimum energy for the spectrum in keV.
    emax : float
        The maximum energy for the spectrum in keV.
    nbins : integer
        The number of channels in the spectrum.
    var_elem : dict
        Abundances of ions. Each dictionary value, specified by the ionic
        symbol, corresponds to the abundance of that symbol. If a float, it is
        understood to be constant and in solar units. If a string or tuple of
        strings, it is assumed to be a spatially varying field.
    binscale : string, optional
        The scale of the energy binning: "linear" or "log".
        Default: "linear"
    temperature_field : string or (ftype, fname) tuple, optional
        The yt temperature field to use for the thermal modeling. Must have
        units of Kelvin. Default: ("gas","temperature")
    emission_measure_field : string or (ftype, fname) tuple, optional
        The yt emission measure field to use for the thermal modeling. Must
        have units of cm^-3. Default: ("gas","emission_measure")
    h_fraction : float, string, or tuple of strings, optional
        The hydrogen mass fraction. If a float, assumes a constant mass
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
    thermal_broad : boolean, optional
        Whether the spectral lines should be thermally
        broadened. Default: True
    model_root : string, optional
        The directory root where the model files are stored. If not provided,
        a default location known to pyXSIM is used.
    model_vers : string, optional
        The version identifier string for the model files, e.g.
        "2.0.2". Default is "3.0.9".
    nolines : boolean, optional
        Turn off lines entirely for generating emission.
        Default: False
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
    prng : integer or :class:`~numpy.random.RandomState` object
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is to use the :mod:`numpy.random` module.

    Examples
    --------
    >>> var_elem = {
    ...     "H^0":  ("flash", "h   "),
    ...     "He^0": ("flash", "he  "),
    ...     "He^1": ("flash", "he1 "),
    ...     "He^2": ("flash", "he2 "),
    ...     "O^0":  ("flash", "o   "),
    ...     "O^1":  ("flash", "o1  "),
    ...     "O^2":  ("flash", "o2  "),
    ...     "O^3":  ("flash", "o3  "),
    ...     "O^4":  ("flash", "o4  "),
    ...     "O^5":  ("flash", "o5  "),
    ...     "O^6":  ("flash", "o6  "),
    ...     "O^7":  ("flash", "o7  "),
    ...     "O^8":  ("flash", "o8  ")
    ... }
    >>> source_model = NEISourceModel(0.1, 10.0, 10000, var_elem)
    """

    _nei = True

    def __init__(
        self,
        emin,
        emax,
        nbins,
        var_elem,
        binscale="linear",
        temperature_field=("gas", "temperature"),
        emission_measure_field=("gas", "emission_measure"),
        h_fraction=None,
        kT_min=0.025,
        kT_max=64.0,
        max_density=None,
        min_entropy=None,
        method="invert_cdf",
        thermal_broad=True,
        model_root=None,
        model_vers=None,
        nolines=False,
        abund_table="angr",
        prng=None,
    ):
        super().__init__(
            "apec",
            emin,
            emax,
            nbins,
            0.0,
            binscale=binscale,
            temperature_field=temperature_field,
            emission_measure_field=emission_measure_field,
            h_fraction=h_fraction,
            kT_min=kT_min,
            kT_max=kT_max,
            max_density=max_density,
            min_entropy=min_entropy,
            var_elem=var_elem,
            method=method,
            thermal_broad=thermal_broad,
            model_root=model_root,
            model_vers=model_vers,
            nolines=nolines,
            abund_table=abund_table,
            prng=prng,
        )

    def _prep_repr(self):
        class_name, strs = super()._prep_repr()
        strs.pop("model")
        strs.pop("Zmet")
        return class_name, strs
