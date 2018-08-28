"""
Photon emission and absoprtion models.
"""
import numpy as np
import h5py

from soxs.spectra import ApecGenerator, \
    get_wabs_absorb, get_tbabs_absorb
from soxs.utils import parse_prng
from pyxsim.utils import mylog
from yt.funcs import get_pbar
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.physical_constants import hcgs, clight

hc = (hcgs*clight).in_units("keV*angstrom").v
# NOTE: XSPEC has hc = 12.39854 keV*A, so there may be slight differences in
# placement of spectral lines due to the above
cl = clight.v

class TableApecModel(ApecGenerator):
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
    >>> apec_model = TableApecModel(0.05, 50.0, 1000, apec_vers="3.0.3",
    ...                             thermal_broad=False)
    """
    def __init__(self, emin, emax, nchan, var_elem=None,
                 model_root=None, model_vers=None, 
                 thermal_broad=True, nolines=False,
                 abund_table="angr", nei=False):
        super(TableApecModel, self).__init__(emin, emax, nchan, var_elem=var_elem,
                                             apec_root=model_root, apec_vers=model_vers, 
                                             broadening=thermal_broad, nolines=nolines,
                                             abund_table=abund_table, nei=nei)
        self.nchan = self.nbins

    def prepare_spectrum(self, zobs):
        """
        Prepare the thermal model for execution given a redshift *zobs* for the spectrum.
        """
        cosmic_spec, metal_spec, var_spec = \
            self._get_table(list(range(self.nT)), zobs, 0.0)
        self.cosmic_spec = YTArray(cosmic_spec, "cm**3/s")
        self.metal_spec = YTArray(metal_spec, "cm**3/s")
        if var_spec is None:
            self.var_spec = var_spec
        else:
            self.var_spec = YTArray(var_spec, "cm**3/s")

    def get_spectrum(self, kT):
        """
        Get the thermal emission spectrum given a temperature *kT* in keV. 
        """
        tindex = np.searchsorted(self.Tvals, kT)-1
        if tindex >= self.Tvals.shape[0]-1 or tindex < 0:
            return (YTArray(np.zeros(self.nchan), "cm**3/s"),)*2
        dT = (kT-self.Tvals[tindex])/self.dTvals[tindex]
        cspec_l = self.cosmic_spec[tindex, :]
        mspec_l = self.metal_spec[tindex, :]
        cspec_r = self.cosmic_spec[tindex+1, :]
        mspec_r = self.metal_spec[tindex+1, :]
        cosmic_spec = cspec_l*(1.-dT)+cspec_r*dT
        metal_spec = mspec_l*(1.-dT)+mspec_r*dT
        var_spec = None
        if self.var_spec is not None:
            vspec_l = self.var_spec[:, tindex, :]
            vspec_r = self.var_spec[:, tindex+1, :]
            var_spec = vspec_l*(1.-dT) + vspec_r*dT
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
        spec = super(TableApecModel, self).get_spectrum(temperature, metallicity, 
                                                        redshift, norm, velocity=velocity,
                                                        elem_abund=elem_abund)
        return YTArray(spec.flux*spec.de, "photons/s/cm**2")

thermal_models = {"apec": TableApecModel}

class AbsorptionModel(object):
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
        detected = np.zeros(n_events, dtype='bool')
        nchunk = n_events // 100
        if nchunk == 0:
            nchunk = n_events
        k = 0
        pbar = get_pbar("Absorbing photons", n_events)
        while k < n_events:
            absorb = self.get_absorb(eobs[k:k+nchunk])
            nabs = absorb.size
            randvec = prng.uniform(size=nabs)
            detected[k:k+nabs] = randvec < absorb
            k += nabs
            pbar.update(k)
        pbar.finish()
        return detected

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
