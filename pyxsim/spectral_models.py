"""
Photon emission and absoprtion models.
"""
import numpy as np
from scipy.interpolate import interp1d
from soxs.constants import K_per_keV
from soxs.spectra import get_tbabs_absorb, get_wabs_absorb
from soxs.thermal_spectra import (
    CIEGenerator,
    CloudyCIEGenerator,
    IGMGenerator,
    MekalGenerator,
)
from soxs.utils import parse_prng, regrid_spectrum
from yt.units.yt_array import YTArray, YTQuantity

from pyxsim.lib.interpolate import interp1d_spec, interp2d_spec


class SpectralInterpolator1D:
    def __init__(self, tbins, cosmic_spec, metal_spec, var_spec):
        self.tbins = tbins.astype("float64")
        self.cosmic_spec = cosmic_spec
        self.metal_spec = metal_spec
        if var_spec is None:
            self.var_spec = np.zeros((1, 1, 1))
            self.do_var = False
        else:
            self.var_spec = var_spec
            self.do_var = True

    def __call__(self, t_vals):
        x_i = (np.digitize(t_vals, self.tbins) - 1).astype("int32")
        if np.any((x_i == -1) | (x_i == len(self.tbins) - 1)):
            x_i = np.minimum(np.maximum(x_i, 0), len(self.tbins) - 2)
        c_vals, m_vals, v_vals = interp1d_spec(
            self.cosmic_spec,
            self.metal_spec,
            self.var_spec,
            t_vals,
            self.tbins,
            x_i,
            self.do_var,
        )
        return c_vals, m_vals, v_vals


class SpectralInterpolator2D:
    def __init__(self, tbins, dbins, cosmic_spec, metal_spec, var_spec):
        self.tbins = tbins.astype("float64")
        self.dbins = dbins.astype("float64")
        self.cosmic_spec = cosmic_spec
        self.metal_spec = metal_spec
        if var_spec is None:
            self.var_spec = np.zeros((1, 1, 1))
            self.do_var = False
        else:
            self.var_spec = var_spec
            self.do_var = True

    def __call__(self, t_vals, d_vals):
        x_i = (np.digitize(t_vals, self.tbins) - 1).astype("int32")
        if np.any((x_i == -1) | (x_i == len(self.tbins) - 1)):
            x_i = np.minimum(np.maximum(x_i, 0), len(self.tbins) - 2)
        y_i = (np.digitize(d_vals, self.dbins) - 1).astype("int32")
        if np.any((y_i == -1) | (y_i == len(self.dbins) - 1)):
            y_i = np.minimum(np.maximum(y_i, 0), len(self.dbins) - 2)
        c_vals, m_vals, v_vals = interp2d_spec(
            self.cosmic_spec,
            self.metal_spec,
            self.var_spec,
            t_vals,
            self.tbins,
            x_i,
            d_vals,
            self.dbins,
            y_i,
            self.do_var,
        )
        return c_vals, m_vals, v_vals


class ThermalSpectralModel:
    _logT = False

    def _Tconv(self, kT):
        if self._logT:
            return np.log10(kT * K_per_keV)
        else:
            return kT

    def get_spectrum(self, kT):
        """
        Get the thermal emission spectrum given a temperature *kT* in keV.
        """
        kT = np.atleast_1d(self._Tconv(kT))
        return self.si(kT)

    def make_fluxf(self, emin, emax, energy=False):
        eidxs = (self.ebins[:-1] > emin) & (self.ebins[1:] < emax)
        emid = self.emid[eidxs]
        if energy:
            cosmic_flux = (self.cosmic_spec[:, eidxs] * emid).sum(axis=-1)
            metal_flux = (self.metal_spec[:, eidxs] * emid).sum(axis=-1)
        else:
            cosmic_flux = self.cosmic_spec[:, eidxs].sum(axis=-1)
            metal_flux = self.metal_spec[:, eidxs].sum(axis=-1)
        cf = interp1d(
            self.Tvals, cosmic_flux, fill_value=0.0, assume_sorted=True, copy=False
        )
        mf = interp1d(
            self.Tvals, metal_flux, fill_value=0.0, assume_sorted=True, copy=False
        )
        if self.var_spec is not None:
            if energy:
                var_flux = (self.var_spec[:, :, eidxs] * emid).sum(axis=-1)
            else:
                var_flux = self.var_spec[:, :, eidxs].sum(axis=-1)
            vf = interp1d(
                self.Tvals,
                var_flux,
                axis=1,
                fill_value=0.0,
                assume_sorted=True,
                copy=False,
            )
        else:

            def vf(kT):
                pass

        def _fluxf(kT):
            kT = self._Tconv(kT)
            return cf(kT), mf(kT), vf(kT)

        return _fluxf


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
    nbins : integer
        The number of bins in the spectral model. If one
        is thermally broadening lines, it is recommended that
        this value result in an energy resolution per channel
        of roughly 1 eV or smaller.
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
        "3.0.3". Default: 3.0.9
    thermal_broad : boolean, optional
        Whether the spectral lines should be thermally
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
        If True, use the non-equilibrium ionization tables. Only available
        for the "apec" model. Default: False

    Examples
    --------
    >>> apec_model = TableCIEModel("apec", 0.05, 50.0, 1000, model_vers="3.0.3",
    ...                            thermal_broad=False)
    """

    def __init__(
        self,
        model,
        emin,
        emax,
        nbins,
        kT_min,
        kT_max,
        binscale="linear",
        var_elem=None,
        model_root=None,
        model_vers=None,
        thermal_broad=True,
        nolines=False,
        abund_table="angr",
        nei=False,
    ):
        self.cgen = CIEGenerator(
            model,
            emin,
            emax,
            nbins,
            binscale=binscale,
            var_elem=var_elem,
            model_root=model_root,
            model_vers=model_vers,
            broadening=thermal_broad,
            nolines=nolines,
            abund_table=abund_table,
            nei=nei,
        )
        self.nbins = self.cgen.nbins
        self.ebins = self.cgen.ebins
        self.emid = self.cgen.emid
        self.var_elem_names = self.cgen.var_elem_names
        self.var_ion_names = self.cgen.var_ion_names
        self.atable = self.cgen.atable
        self.de = np.diff(self.ebins)
        self.kT_min = kT_min
        self.kT_max = kT_max
        self.idx_min = max(np.searchsorted(self.cgen.Tvals, kT_min) - 1, 0)
        self.idx_max = min(
            np.searchsorted(self.cgen.Tvals, kT_max) + 1, self.cgen.nT - 1
        )
        self.Tvals = self.cgen.Tvals[self.idx_min : self.idx_max]
        self.nT = self.Tvals.size
        self.dTvals = np.diff(self.Tvals)
        self.model_vers = self.cgen.model_vers
        self.model_root = self.cgen.model_root

    def prepare_spectrum(self, zobs):
        """
        Prepare the thermal model for execution given a redshift *zobs* for the spectrum.
        """
        cosmic_spec, metal_spec, var_spec = self.cgen._get_table(
            list(range(self.idx_min, self.idx_max)), zobs, 0.0
        )
        self.cosmic_spec = cosmic_spec
        self.metal_spec = metal_spec
        self.var_spec = var_spec
        self.si = SpectralInterpolator1D(
            self.Tvals, self.cosmic_spec, self.metal_spec, self.var_spec
        )


class Atable1DSpectralModel(ThermalSpectralModel):
    _logT = True

    def __init__(self, sgen):
        self.sgen = sgen
        self.nbins = self.sgen.nbins
        self.ebins = self.sgen.ebins
        self.emid = self.sgen.emid
        self.var_elem = self.sgen.var_elem
        self.var_elem_names = self.sgen.var_elem
        self.atable = self.sgen.atable
        self.de = self.sgen.de
        self.binscale = self.sgen.binscale
        self.Tvals = self.sgen.Tvals

    def prepare_spectrum(self, zobs):
        eidxs, ne, ebins, emid, de = self.sgen._get_energies(zobs)
        cosmic_spec, metal_spec, var_spec = self.sgen._get_table(ne, eidxs, zobs)
        self.cosmic_spec = 1.0e-14 * regrid_spectrum(self.ebins, ebins, cosmic_spec)
        self.metal_spec = 1.0e-14 * regrid_spectrum(self.ebins, ebins, metal_spec)
        if var_spec is not None:
            var_spec = 1.0e-14 * regrid_spectrum(self.ebins, ebins, var_spec)
        self.var_spec = var_spec
        self.si = SpectralInterpolator1D(
            self.Tvals, self.cosmic_spec, self.metal_spec, self.var_spec
        )


class MekalSpectralModel(Atable1DSpectralModel):
    def __init__(
        self, emin, emax, nbins, binscale="linear", var_elem=None, abund_table="angr"
    ):
        mgen = MekalGenerator(
            emin,
            emax,
            nbins,
            binscale=binscale,
            var_elem=var_elem,
            abund_table=abund_table,
        )
        super().__init__(mgen)
        self.var_ion_names = []


class CloudyCIESpectralModel(Atable1DSpectralModel):
    def __init__(
        self,
        emin,
        emax,
        nbins,
        binscale="linear",
        var_elem=None,
        model_vers=None,
    ):
        cgen = CloudyCIEGenerator(
            emin,
            emax,
            nbins,
            binscale=binscale,
            var_elem=var_elem,
            model_vers=model_vers,
        )
        super().__init__(cgen)
        self.var_ion_names = []
        self.model_vers = cgen.model_vers


class IGMSpectralModel(ThermalSpectralModel):
    _logT = True
    r"""
    A spectral model for a thermal plasma including photoionization and
    resonant scattering from the CXB based on Khabibullin & Churazov 2019
    (https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4972K/) and Churazov
    et al. 2001 (https://ui.adsabs.harvard.edu/abs/2001MNRAS.323...93C/).

    For temperatures higher than kT ~ 1.09 keV, a Cloudy-based CIE model
    is used to compute the spectrum.

    Assumes the abundance tables from Feldman 1992.

    Table data and README files can be found at
    https://wwwmpa.mpa-garching.mpg.de/~ildar/igm/v2x/.

    Parameters
    ----------
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The minimum energy for the spectral model.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The maximum energy for the spectral model.
    nbins : integer
        The number of bins in the spectral model.
    resonant_scattering : boolean, optional
        Whether or not to include the effects of resonant scattering
        from CXB photons. Default: False
    cxb_factor : float, optional
        The fraction of the CXB photons that are resonant scattered to enhance
        the lines. Default: 0.5
    var_elem : list of strings, optional
        The names of elements to allow to vary freely
        from the single abundance parameter. Default:
        None
    """

    def __init__(
        self,
        emin,
        emax,
        nbins,
        binscale="linear",
        resonant_scattering=False,
        cxb_factor=0.5,
        var_elem=None,
        model_vers=None,
    ):
        self.igen = IGMGenerator(
            emin,
            emax,
            nbins,
            binscale=binscale,
            resonant_scattering=resonant_scattering,
            cxb_factor=cxb_factor,
            var_elem=var_elem,
            model_vers=model_vers,
        )
        self.ebins = self.igen.ebins
        self.emid = self.igen.emid
        self.de = self.igen.de
        self.nbins = nbins
        self.var_elem = self.igen.var_elem
        self.nvar_elem = self.igen.nvar_elem
        self.min_table_kT = 10 ** self.igen.Tvals[0] / K_per_keV
        self.max_table_kT = 10 ** self.igen.Tvals[-1] / K_per_keV
        self.min_table_nH = 10 ** self.igen.Dvals[0]
        self.max_table_nH = 10 ** self.igen.Dvals[-1]
        self.Tvals = self.igen.Tvals
        self.Dvals = self.igen.Dvals
        self.dTvals = self.igen.dTvals
        self.dDvals = self.igen.dDvals
        self.n_T = self.igen.n_T
        self.n_D = self.igen.n_D
        self.binscale = self.igen.binscale
        self.cie_model = CloudyCIESpectralModel(
            emin,
            emax,
            nbins,
            binscale=self.binscale,
            var_elem=self.var_elem,
            model_vers=model_vers,
        )
        self.model_vers = self.igen.model_vers

    def prepare_spectrum(self, zobs):
        """
        Prepare the thermal model for execution given a redshift *zobs* for the spectrum.
        """
        eidxs, ne, ebins, emid, de = self.igen._get_energies(zobs)
        cosmic_spec, metal_spec, var_spec = self.igen._get_table(ne, eidxs, zobs)
        self.cosmic_spec = 1.0e-14 * regrid_spectrum(self.ebins, ebins, cosmic_spec)
        self.metal_spec = 1.0e-14 * regrid_spectrum(self.ebins, ebins, metal_spec)
        if var_spec is not None:
            var_spec = 1.0e-14 * regrid_spectrum(self.ebins, ebins, var_spec)
        self.var_spec = var_spec
        self.cie_model.prepare_spectrum(zobs)
        self.si = SpectralInterpolator2D(
            self.Tvals, self.Dvals, self.cosmic_spec, self.metal_spec, self.var_spec
        )

    def get_spectrum(self, kT, nH):
        kT = np.atleast_1d(kT)
        nH = np.atleast_1d(nH)
        use_igm = (kT >= self.min_table_kT) & (kT <= self.max_table_kT)
        use_igm &= (nH >= self.min_table_nH) & (nH <= self.max_table_nH)
        use_cie = ~use_igm
        cspec = np.zeros((kT.size, self.nbins))
        mspec = np.zeros((kT.size, self.nbins))
        if self.var_spec is not None:
            vspec = np.zeros((self.nvar_elem, kT.size, self.nbins))
        else:
            vspec = None
        n_cie = use_cie.sum()
        n_igm = use_igm.sum()
        if n_igm > 0:
            nHi = nH[use_igm]
            c1, m1, v1 = self._get_spectrum_2d(kT[use_igm], nHi)
            cspec[use_igm, :] = c1 / nHi[:, np.newaxis]
            mspec[use_igm, :] = m1 / nHi[:, np.newaxis]
            if self.var_spec is not None:
                vspec[:, use_igm, :] = v1 / nHi[np.newaxis, :, np.newaxis]
        if n_cie > 0:
            c2, m2, v2 = self.cie_model.get_spectrum(kT[use_cie])
            cspec[use_cie, :] = c2
            mspec[use_cie, :] = m2
            if self.var_spec is not None:
                vspec[:, use_cie, :] = v2
        return cspec, mspec, vspec

    def _get_spectrum_2d(self, kT, nH):
        lkT = np.atleast_1d(np.log10(kT * K_per_keV))
        lnH = np.atleast_1d(np.log10(nH))
        cspec, mspec, vspec = self.si(lkT, lnH)
        return cspec, mspec, vspec

    def make_fluxf(self, emin, emax, energy=False):
        eidxs = (self.ebins[:-1] > emin) & (self.ebins[1:] < emax)
        emid = self.emid[eidxs]
        if energy:
            cosmic_flux = (self.cosmic_spec[:, eidxs] * emid).sum(axis=-1)
            metal_flux = (self.metal_spec[:, eidxs] * emid).sum(axis=-1)
        else:
            cosmic_flux = self.cosmic_spec[:, eidxs].sum(axis=-1)
            metal_flux = self.metal_spec[:, eidxs].sum(axis=-1)
        if self.var_spec is not None:
            if energy:
                var_flux = (self.var_spec[:, :, eidxs] * emid).sum(axis=-1)
            else:
                var_flux = self.var_spec[:, :, eidxs].sum(axis=-1)
        else:
            var_flux = None
        cie_fluxf = self.cie_model.make_fluxf(emin, emax, energy=energy)

        def _fluxf(kT, nH):
            kT = np.atleast_1d(kT)
            nH = np.atleast_1d(nH)
            use_igm = (kT >= self.min_table_kT) & (kT <= self.max_table_kT)
            use_igm &= (nH >= self.min_table_nH) & (nH <= self.max_table_nH)
            use_cie = ~use_igm
            cflux = np.zeros(kT.size)
            mflux = np.zeros(kT.size)
            if self.var_spec is not None:
                vflux = np.zeros((self.nvar_elem, kT.size))
            else:
                vflux = None
            n_cie = use_cie.sum()
            n_igm = use_igm.sum()
            if n_igm > 0:
                nHi = nH[use_igm]
                c1, m1, v1 = self._get_flux_2d(
                    kT[use_igm], nHi, cosmic_flux, metal_flux, var_flux
                )
                cflux[use_igm] = c1 / nHi
                mflux[use_igm] = m1 / nHi
                if self.var_spec is not None:
                    vflux[:, use_igm] = v1 / nHi[np.newaxis, :]
            if n_cie > 0:
                c2, m2, v2 = cie_fluxf(kT[use_cie])
                cflux[use_cie] = c2
                mflux[use_cie] = m2
                if self.var_spec is not None:
                    vflux[:, use_cie] = v2
            return cflux, mflux, vflux

        return _fluxf

    def _get_flux_2d(self, kT, nH, cf, mf, vf):
        lkT = np.atleast_1d(np.log10(kT * K_per_keV))
        lnH = np.atleast_1d(np.log10(nH))
        tidxs = np.searchsorted(self.Tvals, lkT) - 1
        didxs = np.searchsorted(self.Dvals, lnH) - 1
        dT = (lkT - self.Tvals[tidxs]) / self.dTvals[tidxs]
        dn = (lnH - self.Dvals[didxs]) / self.dDvals[didxs]
        idx1 = np.ravel_multi_index((didxs + 1, tidxs + 1), (self.n_D, self.n_T))
        idx2 = np.ravel_multi_index((didxs + 1, tidxs), (self.n_D, self.n_T))
        idx3 = np.ravel_multi_index((didxs, tidxs + 1), (self.n_D, self.n_T))
        idx4 = np.ravel_multi_index((didxs, tidxs), (self.n_D, self.n_T))
        dx1 = dT * dn
        dx2 = dn - dx1
        dx3 = dT - dx1
        dx4 = 1.0 + dx1 - dT - dn
        cflux = dx1 * cf[idx1] + dx2 * cf[idx2] + dx3 * cf[idx3] + dx4 * cf[idx4]
        mflux = dx1 * mf[idx1] + dx2 * mf[idx2] + dx3 * mf[idx3] + dx4 * mf[idx4]
        if vf is not None:
            vflux = dx1[np.newaxis, :] * vf[:, idx1]
            vflux += dx2[np.newaxis, :] * vf[:, idx2]
            vflux += dx3[np.newaxis, :] * vf[:, idx3]
            vflux += dx4[np.newaxis, :] * vf[:, idx4]
        else:
            vflux = None
        return cflux, mflux, vflux


class AbsorptionModel:
    _name = ""

    def __init__(self, nH, energy, cross_section):
        self.nH = YTQuantity(nH * 1.0e22, "cm**-2")
        self.emid = YTArray(energy, "keV")
        self.sigma = YTArray(cross_section, "cm**2")

    def get_absorb(self, e):
        """
        Get the absorption spectrum.
        """
        sigma = np.interp(e, self.emid, self.sigma, left=0.0, right=0.0)
        return np.exp(-sigma * self.nH)

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
            return np.array([], dtype="bool")
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

    def __init__(self, nH, abund_table="angr"):
        self.nH = YTQuantity(nH, "1.0e22*cm**-2")
        self.abund_table = abund_table

    def get_absorb(self, e):
        e = np.array(e)
        return get_tbabs_absorb(e, self.nH.v, abund_table=self.abund_table)


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

    def __init__(self, nH, abund_table="angr"):
        self.nH = YTQuantity(nH, "1.0e22*cm**-2")
        self.abund_table = abund_table

    def get_absorb(self, e):
        e = np.array(e)
        return get_wabs_absorb(e, self.nH.v)


absorb_models = {"wabs": WabsModel, "tbabs": TBabsModel}
