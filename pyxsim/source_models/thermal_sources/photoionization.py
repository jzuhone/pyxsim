from pyxsim.spectral_models import PionSpectralModel

from .base import ThermalSourceModel


class PionSourceModel(ThermalSourceModel):
    """
    A source model for a thermal plasma including photoionization and
    resonant scattering from the CXB based on Khabibullin & Churazov 2019
    (https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4972K/) and Churazov
    et al. 2001 (https://ui.adsabs.harvard.edu/abs/2001MNRAS.323...93C/).

    For temperatures higher than kT ~ 1.09 keV, a Cloudy-based CIE model
    is used to compute the spectrum.

    Assumes the abundance table from Feldman 1992.

    Table data and README files can be found at
    https://wwwmpa.mpa-garching.mpg.de/~ildar/igm/v3/.

    Parameters
    ----------
    emin : float, (value, unit) tuple, unyt_quantity, or Quantity
        The minimum energy for the spectral model.
    emax : float, (value, unit) tuple, unyt_quantity, or Quantity
        The maximum energy for the spectral model.
    nbins : integer
        The number of bins in the spectral model. If one
        is thermally broadening lines, it is recommended that
        this value result in an energy resolution per channel
        of roughly 1 eV or smaller.
    Zmet : float, string, or tuple of strings
        The metallicity. If a float, assumes a constant metallicity throughout
        in solar units. If a string or tuple of strings, is taken to be the
        name of the metallicity field.
    binscale : string, optional
        The scale of the energy binning: "linear" or "log".
        Default: "linear"
    resonant_scattering : boolean, optional
        Whether to include the effects of resonant scattering
        from CXB photons. Default: False
    cxb_factor : float, optional
        The fraction of the CXB photons that are resonant scattered to enhance
        the lines. Default: 0.5
    model_vers : string, optional
        The version identifier string for the model files, "3" or "4".
    nh_field : string or (ftype, fname) tuple, optional
        The yt hydrogen nuclei density field (meaning all hydrogen, ionized or not)
        to use for the model. Must have units of cm**-3.
        Default: ("gas", "H_nuclei_density")
    temperature_field : string or (ftype, fname) tuple, optional
        The yt temperature field to use for the thermal modeling. Must have
        units of Kelvin. Default: ("gas", "temperature")
    emission_measure_field : string or (ftype, fname) tuple, optional
        The yt emission measure field to use for the thermal modeling. Must
        have units of cm^-3. Default: ("gas", "emission_measure")
    h_fraction : float, string, or tuple of strings, optional
        The hydrogen mass fraction. If a float, assumes a constant mass
        fraction of hydrogen throughout. If a string or tuple of strings,
        is taken to be the name of the hydrogen fraction field. Defaults to
        the appropriate value for the Feldman abundance table.
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
        assumed to be a spatially varying field. Default: None
    method : string, optional
        The method used to generate the photon energies from the spectrum:
        "invert_cdf": Invert the cumulative distribution function of the spectrum.
        "accept_reject": Acceptance-rejection method using the spectrum.
        The first method should be sufficient for most cases.
    model_vers : string, optional
        The version of the Pion tables to use in the calculations.
        Options are:
        "4_lo": Tables computed from Cloudy using a continuum resolution
        of 0.1 with a range of 0.05 to 10 keV.
        "4_hi": Tables computed from Cloudy using enhanced continuum
        resolution of 0.025 with a range of 0.05 to 10 keV. Excellent
        energy resolution, but may be expensive to evaluate.
        Default: "4_lo"
    prng : integer or :class:`~numpy.random.RandomState` object
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is to use the :mod:`numpy.random` module.
    """

    _nei = False
    _density_dependence = True

    def __init__(
        self,
        emin,
        emax,
        nbins,
        Zmet,
        binscale="linear",
        resonant_scattering=False,
        cxb_factor=0.5,
        nh_field=("gas", "H_nuclei_density"),
        temperature_field=("gas", "temperature"),
        emission_measure_field=("gas", "emission_measure"),
        h_fraction=None,
        kT_min=0.025,
        kT_max=64.0,
        max_density=None,
        min_entropy=None,
        var_elem=None,
        method="invert_cdf",
        model_vers="4_lo",
        prng=None,
    ):
        var_elem_keys = list(var_elem.keys()) if var_elem else None
        spectral_model = PionSpectralModel(
            emin,
            emax,
            nbins,
            binscale=binscale,
            resonant_scattering=resonant_scattering,
            cxb_factor=cxb_factor,
            var_elem=var_elem_keys,
            model_vers=model_vers,
        )
        nH_min = 10 ** spectral_model.Dvals[0]
        nH_max = 10 ** spectral_model.Dvals[-1]
        super().__init__(
            spectral_model,
            emin,
            emax,
            nbins,
            Zmet,
            binscale=binscale,
            kT_min=kT_min,
            kT_max=kT_max,
            nH_min=nH_min,
            nH_max=nH_max,
            var_elem=var_elem,
            max_density=max_density,
            min_entropy=min_entropy,
            method=method,
            abund_table="feld",
            prng=prng,
            temperature_field=temperature_field,
            h_fraction=h_fraction,
            emission_measure_field=emission_measure_field,
        )
        self.nh_field = nh_field
        self.resonant_scattering = resonant_scattering
        self.cxb_factor = cxb_factor

    def _prep_repr(self):
        class_name, strs = super()._prep_repr()
        strs["nh_field"] = self.nh_field
        strs["resonant_scattering"] = self.resonant_scattering
        strs["cxb_factor"] = self.cxb_factor
        return class_name, strs


class IGMSourceModel(PionSourceModel):
    pass
