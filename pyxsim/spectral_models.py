"""
Photon emission and absoprtion models.
"""
import numpy as np

from soxs.thermal_spectra import CIEGenerator, IGMGenerator
from soxs.spectra import \
    get_wabs_absorb, get_tbabs_absorb
from soxs.constants import K_per_keV
from soxs.utils import parse_prng
from yt.units.yt_array import YTArray, YTQuantity
from scipy.interpolate import interp1d


class ThermalSpectralModel:
    pass


class TableCIEModel(ThermalSpectralModel):
    r"""
    Initialize a thermal gas emission model from the AtomDB APEC tables
    available at http://www.atomdb.org. This code borrows heavily from Python
    routines used to read the APEC tables developed by Adam Foster at the
    CfA (afoster@cfa.harvard.edu).

    Parameters
    ----------
    emin : float
        The minimum energy for the spectral model.
    emax : float
        The maximum energy for the spectral model.
    nchan : integer
        The number of channels in the spectral model. If one
        is thermally broadening lines, it is recommended that 
        this value result in an energy resolution per channel
        of roughly 1 eV.
    binscale : string, optional
        The scale of the energy binning: "linear" or "log". 
        Default: "linear"
    var_elem : list of strings, optional
        The names of elements to allow to vary freely
        from the single abundance parameter. Default:
        None
    model_root : string, optional
        The directory root where the model files are stored. If 
        not provided, the default SOXS-provided files are used.
    model_vers : string, optional
        The version identifier string for the APEC files, e.g.
        "3.0.3". Default: 3.0.8
    thermal_broad : boolean, optional
        Whether or not the spectral lines should be thermally
        broadened. Default: True
    nolines : boolean, optional
        Turn off lines entirely for generating spectra.
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
    nei : boolean, optional
        If True, use the non-equilibrium ionization tables. These are
        not supplied with pyXSIM/SOXS but must be downloaded separately, in
        which case the *apec_root* parameter must also be set to their
        location. Default: False

    Examples
    --------
    >>> apec_model = TableCIEModel("apec", 0.05, 50.0, 1000, model_vers="3.0.3",
    ...                            thermal_broad=False)
    """
    def __init__(self, model, emin, emax, nchan, binscale="linear", var_elem=None,
                 model_root=None, model_vers=None, 
                 thermal_broad=True, nolines=False,
                 abund_table="angr", nei=False):
        self.cgen = CIEGenerator(model, emin, emax, nchan, binscale=binscale, 
                                 var_elem=var_elem, model_root=model_root, 
                                 model_vers=model_vers, broadening=thermal_broad, 
                                 nolines=nolines, abund_table=abund_table, nei=nei)
        self.nchan = self.cgen.nbins
        self.ebins = self.cgen.ebins
        self.emid = self.cgen.emid
        self.var_elem_names = self.cgen.var_elem_names
        self.var_ion_names = self.cgen.var_ion_names
        self.atable = self.cgen.atable
        self.de = np.diff(self.ebins)

    def prepare_spectrum(self, zobs, kT_min, kT_max):
        """
        Prepare the thermal model for execution given a redshift *zobs* for the spectrum.
        """
        idx_min = max(np.searchsorted(self.cgen.Tvals, kT_min)-1, 0)
        idx_max = min(np.searchsorted(self.cgen.Tvals, kT_max)+1, self.cgen.nT-1)
        cosmic_spec, metal_spec, var_spec = \
            self.cgen._get_table(list(range(idx_min, idx_max)), zobs, 0.0)
        self.Tvals = self.cgen.Tvals[idx_min:idx_max]
        self.nT = self.Tvals.size
        self.dTvals = np.diff(self.Tvals)
        self.cosmic_spec = cosmic_spec
        self.metal_spec = metal_spec
        self.var_spec = var_spec
        self.cf = interp1d(self.Tvals, self.cosmic_spec, axis=0, fill_value=0.0,
                           assume_sorted=True, copy=False)
        self.mf = interp1d(self.Tvals, self.metal_spec, axis=0, fill_value=0.0,
                           assume_sorted=True, copy=False)
        if var_spec is not None:
            self.vf = interp1d(self.Tvals, self.var_spec, axis=1, fill_value=0.0,
                               assume_sorted=True, copy=False)
        else:
            self.vf = None

    def get_spectrum(self, kT):
        """
        Get the thermal emission spectrum given a temperature *kT* in keV. 
        """
        kT = np.atleast_1d(kT)
        var_spec = None
        cosmic_spec = self.cf(kT)
        metal_spec = self.mf(kT)
        if self.var_spec is not None:
            var_spec = self.vf(kT)
        return cosmic_spec, metal_spec, var_spec

    def return_spectrum(self, temperature, metallicity, redshift, norm,
                        velocity=0.0, elem_abund=None):
        """
        Given the properties of a thermal plasma, return a spectrum.

        Parameters
        ----------
        temperature : float
            The temperature of the plasma in keV.
        metallicity : float
            The metallicity of the plasma in solar units.
        redshift : float
            The redshift of the plasma.
        norm : float
            The normalization of the model, in the standard Xspec units of
            1.0e-14*EM/(4*pi*(1+z)**2*D_A**2).
        velocity : float, optional
            Velocity broadening parameter in km/s. Default: 0.0
        elem_abund : dict of element name, float pairs
            A dictionary of elemental abundances to vary
            freely of the abund parameter. Default: None
        """
        spec = self.cgen.get_spectrum(temperature, metallicity, redshift, norm, 
                                      velocity=velocity, elem_abund=elem_abund)
        return YTArray(spec.flux*spec.de, "photons/s/cm**2")


class IGMSpectralModel(ThermalSpectralModel):
    r"""
    A spectral model for a thermal plasma including photoionization and 
    resonant scattering from the CXB based on Khabibullin & Churazov 2019
    (https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4972K/) and Churazov 
    et al. 2001 (https://ui.adsabs.harvard.edu/abs/2001MNRAS.323...93C/).

    For temperatures higher than kT ~ 1.09 keV, APEC is used to compute the
    spectrum. 
 
    Assumes the abundance tables from Feldman 1992.

    Table data and README files can be found at
    https://wwwmpa.mpa-garching.mpg.de/~ildar/igm/v2x/.

    Parameters
    ----------
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The minimum energy for the spectral model.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The maximum energy for the spectral model.
    resonant_scattering : boolean, optional
        Whether or not to include the effects of resonant scattering
        from CXB photons. Default: False
    cxb_factor : float, optional
        The fraction of the CXB photons that are resonant scattered to enhance
        the lines. Default: 0.5
    var_elem_option: integer, optional
        An integer to choose between options for variable elements, which are:
        1: specify abundances of O, Ne, and Fe separately from other metals
        2: specify abundances of O, Ne, Mg, Si, S, and Fe separately from other
           metals
        Default: None, which means no metal abundances can be specified
        separately.
    var_elem : list of strings, optional
        The names of elements to allow to vary freely
        from the single abundance parameter. Default:
        None
    """
    def __init__(self, emin, emax, resonant_scattering=False, cxb_factor=0.5,
                 var_elem_option=None, var_elem=None):
        self.igen = IGMGenerator(emin, emax, resonant_scattering=resonant_scattering,
                                 cxb_factor=cxb_factor, var_elem_option=var_elem_option)
        if var_elem is not None:
            if set(var_elem) != set(self.igen.var_elem):
                raise RuntimeError("The supplied set of abundances does not match "
                                   "what is available for 'var_elem_option = "
                                   f"{self.igen.var_elem_option}!\n"
                                   "Free elements: %s\nAbundances: %s" % (set(var_elem),
                                                                          set(self.igen.var_elem)))

        self.var_elem = var_elem
        self.nvar_elem = self.igen.nvar_elem
        self.min_table_kT = 10**self.igen.Tvals[0] / K_per_keV
        self.max_table_kT = 10**self.igen.Tvals[-1] / K_per_keV
        self.min_table_nH = 10**self.igen.Dvals[0]
        self.max_table_nH = 10**self.igen.Dvals[-1]
        self.Tvals = self.igen.Tvals
        self.Dvals = self.igen.Dvals
        self.dTvals = self.igen.dTvals
        self.dDvals = self.igen.dDvals
        self.n_T = self.igen.n_T
        self.n_D = self.igen.n_D

    def prepare_spectrum(self, zobs, kT_min, kT_max):
        """
        Prepare the thermal model for execution given a redshift *zobs* for the spectrum.
        """
        eidxs, self.ne, self.ebins, self.emid, self.de = self.igen._get_energies(zobs)
        self.apec_model = TableCIEModel("apec", self.ebins[0], self.ebins[-1], self.ne,
                                        binscale="log", var_elem=self.var_elem,
                                        abund_table="feld")
        cosmic_spec, metal_spec, var_spec = self.igen._get_table(self.ne, eidxs, zobs)
        self.cosmic_spec = 1.0e-14*cosmic_spec
        self.metal_spec = 1.0e-14*metal_spec
        if var_spec is not None:
            var_spec *= 1.0e-14
        self.var_spec = var_spec
        self.apec_model.prepare_spectrum(zobs, kT_min, kT_max)

    def get_spectrum(self, kT, nH):
        kT = np.atleast_1d(kT)
        nH = np.atleast_1d(nH)
        use_igm = (kT >= self.min_table_kT) & (kT <= self.max_table_kT)
        use_igm &= (nH >= self.min_table_nH) & (nH <= self.max_table_nH)
        use_apec = ~use_igm
        cspec = np.zeros((kT.size, self.ne))
        mspec = np.zeros((kT.size, self.ne))
        if self.var_spec is not None:
            vspec = np.zeros((self.nvar_elem, kT.size, self.ne))
        else:
            vspec = None
        n_apec = use_apec.sum()
        n_igm = use_igm.sum()
        if n_igm > 0:
            nHi = nH[use_igm]
            c1, m1, v1 = self._get_spectrum_2d(kT[use_igm], nHi)
            cspec[use_igm, :] = c1/nHi[:,np.newaxis]
            mspec[use_igm, :] = m1/nHi[:,np.newaxis]
            if self.var_spec is not None:
                vspec[:,use_igm,:] = v1/nHi[np.newaxis,:,np.newaxis]
        if n_apec > 0:
            c2, m2, v2 = self.apec_model.get_spectrum(kT[use_apec])
            cspec[use_apec, :] = c2
            mspec[use_apec, :] = m2
            if self.var_spec is not None:
                vspec[:,use_apec,:] = v2
        return cspec, mspec, vspec
        
    def _get_spectrum_2d(self, kT, nH):
        lkT = np.atleast_1d(np.log10(kT*K_per_keV))
        lnH = np.atleast_1d(np.log10(nH))
        tidxs = np.searchsorted(self.Tvals, lkT)-1
        didxs = np.searchsorted(self.Dvals, lnH)-1
        dT = (lkT - self.Tvals[tidxs]) / self.dTvals[tidxs]
        dn = (lnH - self.Dvals[didxs]) / self.dDvals[didxs]
        idx1 = np.ravel_multi_index((didxs+1,tidxs+1), (self.n_D, self.n_T))
        idx2 = np.ravel_multi_index((didxs+1,tidxs), (self.n_D, self.n_T))
        idx3 = np.ravel_multi_index((didxs,tidxs+1), (self.n_D, self.n_T))
        idx4 = np.ravel_multi_index((didxs,tidxs), (self.n_D, self.n_T))
        dx1 = dT*dn
        dx2 = dn-dx1
        dx3 = dT-dx1
        dx4 = 1.0+dx1-dT-dn
        cspec = dx1[:,np.newaxis]*self.cosmic_spec[idx1,:]
        cspec += dx2[:,np.newaxis]*self.cosmic_spec[idx2,:]
        cspec += dx3[:,np.newaxis]*self.cosmic_spec[idx3,:]
        cspec += dx4[:,np.newaxis]*self.cosmic_spec[idx4,:]
        mspec = dx1[:,np.newaxis]*self.metal_spec[idx1,:]
        mspec += dx2[:,np.newaxis]*self.metal_spec[idx2,:]
        mspec += dx3[:,np.newaxis]*self.metal_spec[idx3,:]
        mspec += dx4[:,np.newaxis]*self.metal_spec[idx4,:]
        if self.var_spec is not None:
            vspec = dx1[np.newaxis,:,np.newaxis]*self.var_spec[:,idx1,:]
            vspec += dx2[np.newaxis,:,np.newaxis]*self.var_spec[:,idx2,:]
            vspec += dx3[np.newaxis,:,np.newaxis]*self.var_spec[:,idx3,:]
            vspec += dx4[np.newaxis,:,np.newaxis]*self.var_spec[:,idx4,:]
        else:
            vspec = None
        return cspec, mspec, vspec


class AbsorptionModel:
    _name = ""

    def __init__(self, nH, energy, cross_section):
        self.nH = YTQuantity(nH*1.0e22, "cm**-2")
        self.emid = YTArray(energy, "keV")
        self.sigma = YTArray(cross_section, "cm**2")

    def get_absorb(self, e):
        """
        Get the absorption spectrum.
        """
        sigma = np.interp(e, self.emid, self.sigma, left=0.0, right=0.0)
        return np.exp(-sigma*self.nH)

    def absorb_photons(self, eobs, prng=None):
        r"""
        Determine which photons will be absorbed by foreground
        galactic absorption.

        Parameters
        ----------
        eobs : array_like
            The energies of the photons in keV.
        prng : integer, :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        """
        prng = parse_prng(prng)
        n_events = eobs.size
        if n_events == 0:
            return np.array([], dtype='bool')
        absorb = self.get_absorb(eobs)
        randvec = prng.uniform(size=n_events)
        return randvec < absorb


class TBabsModel(AbsorptionModel):
    r"""
    Initialize a Tuebingen-Boulder (Wilms, J., Allen, A., & 
    McCray, R. 2000, ApJ, 542, 914) ISM absorption model.

    Parameters
    ----------
    nH : float
        The foreground column density in units of 10^22 cm^{-2}.

    Examples
    --------
    >>> tbabs_model = TBabsModel(0.1)
    """
    _name = "tbabs"

    def __init__(self, nH):
        self.nH = YTQuantity(nH, "1.0e22*cm**-2")

    def get_absorb(self, e):
        e = np.array(e)
        return get_tbabs_absorb(e, self.nH.v)


class WabsModel(AbsorptionModel):
    r"""
    Initialize a Wisconsin (Morrison and McCammon; ApJ 270, 119) 
    absorption model.

    Parameters
    ----------
    nH : float
        The foreground column density in units of 10^22 cm^{-2}.

    Examples
    --------
    >>> wabs_model = WabsModel(0.1)
    """
    _name = "wabs"

    def __init__(self, nH):
        self.nH = YTQuantity(nH, "1.0e22*cm**-2")

    def get_absorb(self, e):
        e = np.array(e)
        return get_wabs_absorb(e, self.nH.v)


absorb_models = {"wabs": WabsModel,
                 "tbabs": TBabsModel}
