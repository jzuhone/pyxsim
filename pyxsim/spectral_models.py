"""
Photon emission and absoprtion models.
"""
import numpy as np
import os
import h5py

from pyxsim.utils import mylog, check_file_location
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.physical_constants import hcgs, clight
from yt.utilities.physical_ratios import erg_per_keV, amu_grams
from pyxsim.cutils import broaden_lines
from yt.utilities.on_demand_imports import _astropy

hc = (hcgs*clight).in_units("keV*angstrom").v
# NOTE: XSPEC has hc = 12.39854 keV*A, so there may be slight differences in
# placement of spectral lines due to the above
cl = clight.v

class ThermalSpectralModel(object):

    def __init__(self, emin, emax, nchan):
        self.emin = YTQuantity(emin, "keV")
        self.emax = YTQuantity(emax, "keV")
        self.nchan = nchan
        self.ebins = YTArray(np.linspace(self.emin, self.emax, nchan+1), "keV")
        self.de = np.diff(self.ebins)
        self.emid = 0.5*(self.ebins[1:]+self.ebins[:-1])

    def prepare_spectrum(self, redshift):
        pass

    def cleanup_spectrum(self):
        pass

    def get_spectrum(self, kT):
        pass

class XSpecThermalModel(ThermalSpectralModel):
    r"""
    Initialize a thermal gas emission model from PyXspec.

    Parameters
    ----------
    model_name : string
        The name of the thermal emission model.
    emin : float
        The minimum energy for the spectral model.
    emax : float
        The maximum energy for the spectral model.
    nchan : integer
        The number of channels in the spectral model.
    thermal_broad : boolean, optional
        Whether or not the spectral lines should be thermally
        broadened.
    settings : dictionary, optional
        A dictionary of key, value pairs (must both be strings)
        that can be used to set various options in XSPEC.

    Examples
    --------
    >>> mekal_model = XSpecThermalModel("mekal", 0.05, 50.0, 1000)
    """
    def __init__(self, model_name, emin, emax, nchan,
                 thermal_broad=False, settings=None):
        mylog.warning("XSpecThermalModel is deprecated and will be removed "
                      "in a future release. Use of TableApecModel is suggested.")
        self.model_name = model_name
        self.thermal_broad = thermal_broad
        if settings is None: settings = {}
        self.settings = settings
        super(XSpecThermalModel, self).__init__(emin, emax, nchan)

    def prepare_spectrum(self, zobs):
        """
        Prepare the thermal model for execution given a redshift *zobs* for the spectrum.
        """
        import xspec
        xspec.Xset.chatter = 0
        if self.thermal_broad:
            xspec.Xset.addModelString("APECTHERMAL","yes")
        for k,v in self.settings.items():
            xspec.Xset.addModelString(k,v)
        xspec.AllModels.setEnergies("%f %f %d lin" %
                                    (self.emin.value, self.emax.value, self.nchan))
        self.model = xspec.Model(self.model_name)
        self.thermal_comp = getattr(self.model, self.model_name)
        if self.model_name == "bremss":
            self.norm = 3.02e-15
        else:
            self.norm = 1.0e-14
        self.thermal_comp.norm = 1.0
        self.thermal_comp.Redshift = zobs

    def get_spectrum(self, kT):
        """
        Get the thermal emission spectrum given a temperature *kT* in keV. 
        """
        self.thermal_comp.kT = kT
        self.thermal_comp.Abundanc = 0.0
        cosmic_spec = np.array(self.model.values(0))
        if self.model_name == "bremss":
            metal_spec = np.zeros(self.nchan)
        else:
            self.thermal_comp.Abundanc = 1.0
            metal_spec = np.array(self.model.values(0)) - cosmic_spec
        cosmic_spec *= self.norm
        metal_spec *= self.norm
        return YTArray(cosmic_spec, "cm**3/s"), YTArray(metal_spec, "cm**3/s")

    def cleanup_spectrum(self):
        del self.thermal_comp
        del self.model

class TableApecModel(ThermalSpectralModel):
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
        The number of channels in the spectral model.
    apec_root : string
        The directory root where the APEC model files are stored. If 
        not provided, the default is to look for them in the pyxsim
        "spectral_files" directory.
    apec_vers : string, optional
        The version identifier string for the APEC files, e.g.
        "2.0.2"
    thermal_broad : boolean, optional
        Whether or not the spectral lines should be thermally
        broadened.

    Examples
    --------
    >>> apec_model = TableApecModel(0.05, 50.0, 1000, apec_vers="3.0",
    ...                             thermal_broad=True)
    """
    def __init__(self, emin, emax, nchan, apec_root=None,
                 apec_vers="2.0.2", thermal_broad=False):
        if apec_root is None:
            self.cocofile = check_file_location("apec_v%s_coco.fits" % apec_vers,
                                                "spectral_files")
            self.linefile = check_file_location("apec_v%s_line.fits" % apec_vers,
                                                "spectral_files")
        else:
            self.cocofile = os.path.join(apec_root, "apec_v%s_coco.fits" % apec_vers)
            self.linefile = os.path.join(apec_root, "apec_v%s_line.fits" % apec_vers)
        if not os.path.exists(self.cocofile) or not os.path.exists(self.linefile):
            raise IOError("Cannot find the APEC files!\n %s\n, %s" % (self.cocofile,
                                                                      self.linefile))
        super(TableApecModel, self).__init__(emin, emax, nchan)
        self.wvbins = hc/self.ebins[::-1].d
        # H, He, and trace elements
        self.cosmic_elem = [1,2,3,4,5,9,11,15,17,19,21,22,23,24,25,27,29,30]
        # Non-trace metals
        self.metal_elem = [6,7,8,10,12,13,14,16,18,20,26,28]
        self.thermal_broad = thermal_broad
        self.A = np.array([0.0,1.00794,4.00262,6.941,9.012182,10.811,
                           12.0107,14.0067,15.9994,18.9984,20.1797,
                           22.9898,24.3050,26.9815,28.0855,30.9738,
                           32.0650,35.4530,39.9480,39.0983,40.0780,
                           44.9559,47.8670,50.9415,51.9961,54.9380,
                           55.8450,58.9332,58.6934,63.5460,65.3800])

        try:
            self.line_handle = _astropy.pyfits.open(self.linefile)
        except IOError:
            mylog.error("LINE file %s does not exist" % self.linefile)
            raise IOError("LINE file %s does not exist" % self.linefile)
        try:
            self.coco_handle = _astropy.pyfits.open(self.cocofile)
        except IOError:
            mylog.error("COCO file %s does not exist" % self.cocofile)
            raise IOError("COCO file %s does not exist" % self.cocofile)

        self.Tvals = self.line_handle[1].data.field("kT")
        self.nT = len(self.Tvals)
        self.dTvals = np.diff(self.Tvals)
        self.minlam = self.wvbins.min()
        self.maxlam = self.wvbins.max()

    def prepare_spectrum(self, zobs):
        """
        Prepare the thermal model for execution given a redshift *zobs* for the spectrum.
        """
        sfac = 1.0/(1.+zobs)

        cosmic_spec = np.zeros((self.nT, self.nchan))
        metal_spec = np.zeros((self.nT, self.nchan))

        for ikT, kT in enumerate(self.Tvals):
            line_fields, coco_fields = self._preload_data(ikT)
            # First do H,He, and trace elements
            for elem in self.cosmic_elem:
                cosmic_spec[ikT,:] += self._make_spectrum(kT, elem, line_fields, coco_fields, sfac)
            # Next do the metals
            for elem in self.metal_elem:
                metal_spec[ikT,:] += self._make_spectrum(kT, elem, line_fields, coco_fields, sfac)

        self.cosmic_spec = YTArray(cosmic_spec, "cm**3/s")
        self.metal_spec = YTArray(metal_spec, "cm**3/s")

    def _make_spectrum(self, kT, element, line_fields, coco_fields, scale_factor, velocity=0.0):

        tmpspec = np.zeros(self.nchan)

        i = np.where((line_fields['element'] == element) &
                     (line_fields['lambda'] > self.minlam) &
                     (line_fields['lambda'] < self.maxlam))[0]

        E0 = hc/line_fields['lambda'][i].astype("float64")*scale_factor
        amp = line_fields['epsilon'][i].astype("float64")
        ebins = self.ebins.d
        de = self.de.d
        emid = self.emid.d
        if self.thermal_broad:
            sigma = 2.*kT*erg_per_keV/(self.A[element]*amu_grams)
            if velocity is not None:
                sigma += 2.0*velocity*velocity
            sigma = E0*np.sqrt(sigma)/cl
            vec = broaden_lines(E0, sigma, amp, ebins)
        else:
            vec = np.histogram(E0, ebins, weights=amp)[0]
        tmpspec += vec

        ind = np.where((coco_fields['Z'] == element) &
                       (coco_fields['rmJ'] == 0))[0]
        if len(ind) == 0:
            return tmpspec
        else:
            ind = ind[0]

        n_cont = coco_fields['N_Cont'][ind]
        e_cont = coco_fields['E_Cont'][ind][:n_cont]
        continuum = coco_fields['Continuum'][ind][:n_cont]

        tmpspec += np.interp(emid, e_cont*scale_factor, continuum)*de/scale_factor

        n_pseudo = coco_fields['N_Pseudo'][ind]
        e_pseudo = coco_fields['E_Pseudo'][ind][:n_pseudo]
        pseudo = coco_fields['Pseudo'][ind][:n_pseudo]

        tmpspec += np.interp(emid, e_pseudo*scale_factor, pseudo)*de/scale_factor

        return tmpspec*scale_factor

    def _preload_data(self, index):
        line_data = self.line_handle[index+2].data
        coco_data = self.coco_handle[index+2].data
        line_fields = ('element', 'lambda', 'epsilon')
        coco_fields = ('Z', 'rmJ', 'N_Cont', 'E_Cont', 'Continuum', 'N_Pseudo',
                       'E_Pseudo', 'Pseudo')
        line_fields = {el: line_data.field(el) for el in line_fields}
        coco_fields = {el: coco_data.field(el) for el in coco_fields}
        return line_fields, coco_fields

    def get_spectrum(self, kT):
        """
        Get the thermal emission spectrum given a temperature *kT* in keV. 
        """
        tindex = np.searchsorted(self.Tvals, kT)-1
        if tindex >= self.Tvals.shape[0]-1 or tindex < 0:
            return (YTArray(np.zeros(self.nchan), "cm**3/s"),)*2
        dT = (kT-self.Tvals[tindex])/self.dTvals[tindex]
        cspec_l = self.cosmic_spec[tindex,:]
        mspec_l = self.metal_spec[tindex,:]
        cspec_r = self.cosmic_spec[tindex+1,:]
        mspec_r = self.metal_spec[tindex+1,:]
        cosmic_spec = cspec_l*(1.-dT)+cspec_r*dT
        metal_spec = mspec_l*(1.-dT)+mspec_r*dT
        return cosmic_spec, metal_spec

    def return_spectrum(self, temperature, metallicity, redshift, norm, velocity=0.0):
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
        """
        velocity = YTQuantity(velocity, "km/s").in_cgs().v
        scale_factor = 1.0/(1.+redshift)

        tindex = np.searchsorted(self.Tvals, temperature)-1
        if tindex >= self.Tvals.shape[0]-1 or tindex < 0:
            return YTArray(np.zeros(self.nchan), "photons/s/cm**2")
        dT = (temperature-self.Tvals[tindex])/self.dTvals[tindex]

        cosmic_spec = np.zeros(self.nchan)
        metal_spec = np.zeros(self.nchan)

        fac = [1.0-dT, dT]

        for i, ikT in enumerate([tindex, tindex+1]):
            line_fields, coco_fields = self._preload_data(ikT)
            kT = self.Tvals[ikT]
            # First do H,He, and trace elements
            for elem in self.cosmic_elem:
                cosmic_spec += fac[i]*self._make_spectrum(kT, elem, line_fields, coco_fields,
                                                          scale_factor, velocity=velocity)
            # Next do the metals
            for elem in self.metal_elem:
                metal_spec += fac[i]*self._make_spectrum(kT, elem, line_fields, coco_fields,
                                                         scale_factor, velocity=velocity)

        tspec = (cosmic_spec+metallicity*metal_spec)
        return YTArray(1.0e14*norm*tspec, "photons/s/cm**2")

class AbsorptionModel(object):
    def __init__(self, nH, emid, sigma):
        self.nH = YTQuantity(nH*1.0e22, "cm**-2")
        self.emid = YTArray(emid, "keV")
        self.sigma = YTArray(sigma, "cm**2")

    def prepare_spectrum(self):
        pass

    def get_absorb(self, e):
        """
        Get the absorption spectrum.
        """
        sigma = np.interp(e, self.emid, self.sigma, left=0.0, right=0.0)
        return np.exp(-sigma*self.nH)

    def cleanup_spectrum(self):
        pass

    def absorb_photons(self, eobs, prng=np.random):
        r"""
        Determine which photons will be absorbed by foreground
        galactic absorption.

        Parameters
        ----------
        eobs : array_like
            The energies of the photons in keV.
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        """
        mylog.info("Absorbing.")
        self.prepare_spectrum()
        absorb = self.get_absorb(eobs)
        randvec = prng.uniform(size=eobs.shape)
        detected = randvec < absorb
        self.cleanup_spectrum()
        return detected

class XSpecAbsorbModel(AbsorptionModel):
    r"""
    Initialize an absorption model from PyXspec.

    Parameters
    ----------
    model_name : string
        The name of the absorption model.
    nH : float
        The foreground column density *nH* in units of 10^22 cm^{-2}.
    emin : float, optional
        The minimum energy for the spectral model.
    emax : float, optional
        The maximum energy for the spectral model.
    nchan : integer, optional
        The number of channels in the spectral model.
    settings : dictionary, optional
        A dictionary of key, value pairs (must both be strings)
        that can be used to set various options in XSPEC.

    Examples
    --------
    >>> abs_model = XSpecAbsorbModel("wabs", 0.1)
    """
    def __init__(self, model_name, nH, emin=0.01, emax=50.0,
                 nchan=100000, settings=None):
        mylog.warning("XSpecAbsorbModel is deprecated and will be removed "
                      "in a future release. Use of the other models is "
                      "suggested.")
        self.model_name = model_name
        self.nH = YTQuantity(nH*1.0e22, "cm**-2")
        if settings is None: settings = {}
        self.settings = settings
        self.emin = emin
        self.emax = emax
        self.nchan = nchan
        ebins = np.linspace(emin, emax, nchan+1)
        self.emid = YTArray(0.5*(ebins[1:]+ebins[:-1]), "keV")

    def prepare_spectrum(self):
        """
        Prepare the absorption model for execution.
        """
        import xspec
        xspec.Xset.chatter = 0
        xspec.AllModels.setEnergies("%f %f %d lin" %
                                    (self.emin, self.emax, self.nchan))
        self.model = xspec.Model(self.model_name+"*powerlaw")
        self.model.powerlaw.norm = self.nchan/(self.emax-self.emin)
        self.model.powerlaw.PhoIndex = 0.0
        for k,v in self.settings.items():
            xspec.Xset.addModelString(k,v)
        m = getattr(self.model, self.model_name)
        m.nH = 1.0
        self.sigma = YTArray(-np.log(self.model.values(0))*1.0e-22, "cm**2")

    def cleanup_spectrum(self):
        del self.model


class TableAbsorbModel(AbsorptionModel):
    r"""
    Initialize an absorption model from a table stored in an HDF5 file.

    Parameters
    ----------
    filename : string
        The name of the table file.
    nH : float
        The foreground column density *nH* in units of 10^22 cm^{-2}.

    Examples
    --------
    >>> abs_model = TableAbsorbModel("tbabs_table.h5", 0.1)
    """
    def __init__(self, filename, nH):
        self.filename = check_file_location(filename, "spectral_files")
        f = h5py.File(self.filename,"r")
        emid = YTArray(0.5*(f["energy"][1:]+f["energy"][:-1]), "keV")
        sigma = YTArray(f["cross_section"][:], "cm**2")
        f.close()
        super(TableAbsorbModel, self).__init__(nH, emid, sigma)

class TBabsModel(TableAbsorbModel):
    r"""
    Initialize a Tuebingen-Boulder (Wilms, J., Allen, A., & 
    McCray, R. 2000, ApJ, 542, 914) ISM absorption model.

    Parameters
    ----------
    nH : float
        The foreground column density *nH* in units of 10^22 cm^{-2}.

    Examples
    --------
    >>> tbabs_model = TBabsModel(0.1)
    """
    def __init__(self, nH):
        super(TBabsModel, self).__init__("tbabs_table.h5", nH)


emx = np.array([0.0, 0.1, 0.284, 0.4, 0.532, 0.707, 0.867,
                1.303, 1.840, 2.471, 3.210, 4.038, 7.111, 8.331, 10.0])
c0 = np.array([17.3, 34.6, 78.1, 71.4, 95.5, 308.9, 120.6, 141.3,
               202.7,342.7,352.2,433.9,629.0,701.2])
c1 = np.array([608.1, 267.9, 18.8, 66.8, 145.8, -380.6, 169.3,
               146.8, 104.7, 18.7, 18.7, -2.4, 30.9, 25.2])
c2 = np.array([-2150., -476.1 ,4.3, -51.4, -61.1, 294.0, -47.7,
               -31.5, -17.0, 0.0, 0.0, 0.75, 0.0, 0.0])

class WabsModel(AbsorptionModel):
    r"""
    Initialize a Wisconsin (Morrison and McCammon; ApJ 270, 119) 
    absorption model.

    Parameters
    ----------
    nH : float
        The foreground column density *nH* in units of 10^22 cm^{-2}.

    Examples
    --------
    >>> wabs_model = WabsModel(0.1)
    """
    def __init__(self, nH):
        self.nH = YTQuantity(nH*1.0e22, "cm**-2")

    def get_absorb(self, e):
        e = np.array(e)
        idxs = np.minimum(np.searchsorted(emx, e)-1, 13)
        sigma = (c0[idxs]+c1[idxs]*e+c2[idxs]*e*e)*1.0e-24/e**3
        return np.exp(-sigma*self.nH)
