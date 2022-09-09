"""
Classes for specific source models
"""
import numpy as np
from tqdm.auto import tqdm
from pyxsim.utils import mylog
from pyxsim.source_models.sources import SourceModel
from yt.data_objects.static_output import Dataset
from yt.units.yt_array import YTQuantity
from pyxsim.spectral_models import TableCIEModel, IGMSpectralModel, \
    CloudyCIESpectralModel, MekalSpectralModel
from pyxsim.utils import parse_value, compute_H_abund
from soxs.utils import parse_prng
from soxs.constants import elem_names, atomic_weights, metal_elem, \
    abund_tables
from yt.utilities.exceptions import YTFieldNotFound
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects, communication_system, parallel_capable
from more_itertools import chunked

comm = communication_system.communicators[-1]


class ParallelProgressBar:
    def __init__(self, title):
        self.title = title
        mylog.info(f"Starting '{title}'")

    def update(self, *args, **kwargs):
        return

    def close(self):
        mylog.info(f"Finishing '{self.title}'")


class DummyProgressBar:
    def __init__(self):
        pass

    def update(self, *args, **kwargs):
        return

    def close(self):
        return


class ThermalSourceModel(SourceModel):
    _density_dependence = False
    _nei = False

    def __init__(self, spectral_model, emin, emax, nbins, Zmet, binscale="linear",
                 kT_min=0.025, kT_max=64.0, var_elem=None, max_density=None, 
                 method="invert_cdf", abund_table="angr", prng=None, 
                 temperature_field=None, emission_measure_field=None, 
                 h_fraction=None, nH_min=None, nH_max=None):
        super().__init__(prng=prng)
        self.spectral_model = spectral_model
        self.emin = parse_value(emin, "keV")
        self.emax = parse_value(emax, "keV")
        self.nbins = nbins
        self.Zmet = Zmet
        if var_elem is None:
            var_elem = {}
            var_elem_keys = None
            self.num_var_elem = 0
        else:
            var_elem_keys = list(var_elem.keys())
            self.num_var_elem = len(var_elem_keys)
        self.var_elem = var_elem
        self.var_elem_keys = var_elem_keys
        if max_density is not None:
            if not isinstance(max_density, YTQuantity):
                if isinstance(max_density, tuple):
                    max_density = YTQuantity(max_density[0], max_density[1])
                else:
                    max_density = YTQuantity(max_density, "g/cm**3")
        if emission_measure_field is None:
            emission_measure_field = ('gas', 'emission_measure')
        if temperature_field is None:
            temperature_field = ('gas', 'temperature')
        self.temperature_field = temperature_field
        self.emission_measure_field = emission_measure_field
        self.density_field = None  # Will be determined later
        self.nh_field = None # Will be set by the subclass
        self.max_density = max_density
        self.tot_num_cells = 0  # Will be determined later
        self.ftype = "gas"
        self.binscale = binscale
        self.abund_table = abund_table
        self.method = method
        self.prng = parse_prng(prng)
        self.kT_min = kT_min
        self.kT_max = kT_max
        mylog.info(f"kT_min = {kT_min} keV")
        mylog.info(f"kT_max = {kT_max} keV")
        self.nH_min = nH_min
        self.nH_max = nH_max
        self.spectral_norm = None
        self.redshift = None
        self.pbar = None
        self.Zconvert = 1.0
        self.mconvert = {}
        self.atable = abund_tables[abund_table].copy()
        if h_fraction is None:
            h_fraction = compute_H_abund(abund_table)
        self.h_fraction = h_fraction

    def setup_model(self, data_source, redshift, spectral_norm, elim=None):
        if isinstance(data_source, Dataset):
            ds = data_source
        else:
            ds = data_source.ds
        try:
            ftype = ds._get_field_info(
                self.emission_measure_field).name[0]
        except YTFieldNotFound:
            raise RuntimeError(f"The {self.emission_measure_field} field is not "
                               "found. If you do not have species fields in "
                               "your dataset, you may need to set "
                               "default_species_fields='ionized' in the call "
                               "to yt.load().")
        self.ftype = ftype
        self.redshift = redshift
        if not self._nei and not isinstance(self.Zmet, float):
            Z_units = str(ds._get_field_info(self.Zmet).units)
            if Z_units in ["dimensionless", "", "code_metallicity"]:
                Zsum = (self.atable*atomic_weights)[metal_elem].sum()
                self.Zconvert = atomic_weights[1]/Zsum
            elif Z_units == "Zsun":
                self.Zconvert = 1.0
            else:
                raise RuntimeError(f"I don't understand metallicity "
                                   f"units of {Z_units}!")
        if self.num_var_elem > 0:
            for key, value in self.var_elem.items():
                if not isinstance(value, float):
                    if "^" in key:
                        elem = key.split("^")[0]
                    else:
                        elem = key
                    n_elem = elem_names.index(elem)
                    m_units = str(ds._get_field_info(value).units)
                    if m_units in ["dimensionless", "", "code_metallicity"]:
                        m = self.atable[n_elem]*atomic_weights[n_elem]
                        self.mconvert[key] = atomic_weights[1]/m
                    elif m_units == "Zsun":
                        self.mconvert[key] = 1.0
                    else:
                        raise RuntimeError(f"I don't understand units of "
                                           f"{m_units} for element {key}!")
        self.density_field = (ftype, "density")
        mylog.info(f"Using emission measure field "
                   f"'{self.emission_measure_field}'.")
        mylog.info(f"Using temperature field "
                   f"'{self.temperature_field}'.")
        if self.nh_field is not None:
            mylog.info(f"Using nH field '{self.nh_field}'.")
        self.spectral_model.prepare_spectrum(redshift, self.kT_min,
                                             self.kT_max)
        self.ebins = self.spectral_model.ebins
        self.de = self.spectral_model.de
        self.emid = self.spectral_model.emid
        self.bin_edges = np.log10(self.ebins) if self.binscale == "log" else self.ebins
        self.nbins = self.emid.size
        self.spectral_norm = spectral_norm

        if isinstance(data_source, Dataset):
            self.pbar = DummyProgressBar()
        else:
            citer = data_source.chunks([], "io")
            num_cells = 0
            for chunk in parallel_objects(citer):
                num_cells += chunk[self.temperature_field].size
            self.tot_num_cells = comm.mpi_allreduce(num_cells)
            if parallel_capable:
                self.pbar = ParallelProgressBar("Processing cells/particles ")
            else:
                self.pbar = tqdm(leave=True, total=self.tot_num_cells,
                                 desc="Processing cells/particles ")

    def make_spectrum(self, data_source, redshift=0.0, dist=None, cosmology=None):
        self.setup_model(data_source, redshift, 1.0)
        spec = np.zeros(self.nbins)
        for chunk in data_source.chunks([], "io"):
            spec += self.process_data("spectrum", chunk)
        return self._make_spectrum(data_source.ds, self.ebins, spec,
                                   redshift, dist, cosmology)

    def make_fluxf(self, emin, emax, energy=False):
        return self.spectral_model.make_fluxf(emin, emax, energy=energy)

    def process_data(self, mode, chunk, fluxf=None):

        orig_shape = chunk[self.temperature_field].shape
        if len(orig_shape) == 0:
            orig_ncells = 0
        else:
            orig_ncells = np.prod(orig_shape)
        if orig_ncells == 0:
            if mode == "photons":
                return
            else:
                return np.array([])

        ret = np.zeros(orig_ncells)

        cut = True

        if self.max_density is not None:
            cut &= np.ravel(chunk[self.density_field]) < self.max_density
        kT = np.ravel(
            chunk[self.temperature_field].to_value("keV", "thermal"))
        cut &= (kT >= self.kT_min) & (kT <= self.kT_max)
        if self.nh_field is not None:
            nH = np.ravel(chunk[self.nh_field].d[cut])
        else:
            nH = None

        num_cells = cut.sum()

        if mode == "photons":
            if num_cells == 0:
                self.pbar.update(orig_ncells)
                return
            else:
                self.pbar.update(orig_ncells-num_cells)
        elif num_cells == 0:
            return np.zeros(orig_shape)

        kT = np.ravel(kT[cut])

        cell_nrm = np.ravel(
            chunk[self.emission_measure_field].d[cut]*self.spectral_norm
        )

        if isinstance(self.h_fraction, float):
            X_H = self.h_fraction
        else:
            X_H = np.ravel(chunk[self.h_fraction].d[cut])

        if self._nei:
            metalZ = np.zeros(num_cells)
            elem_keys = self.var_ion_keys
        else:
            elem_keys = self.var_elem_keys
            if isinstance(self.Zmet, float):
                metalZ = self.Zmet*np.ones(num_cells)
            else:
                mZ = chunk[self.Zmet]
                fac = self.Zconvert
                if str(mZ.units) != "Zsun":
                    fac /= X_H
                metalZ = np.ravel(mZ.d[cut]*fac)

        elemZ = None
        if self.num_var_elem > 0:
            elemZ = np.zeros((self.num_var_elem, num_cells))
            for j, key in enumerate(elem_keys):
                value = self.var_elem[key]
                if isinstance(value, float):
                    elemZ[j,:] = value
                else:
                    eZ = chunk[value]
                    fac = self.mconvert[key]
                    if str(eZ.units) != "Zsun":
                        fac /= X_H
                    elemZ[j,:] = np.ravel(eZ.d[cut]*fac)

        if self.observer == "internal" and mode == "photons":
            pos = np.array([np.ravel(chunk[self.p_fields[i]].to_value("kpc")[cut])
                            for i in range(3)])
            r2 = self.compute_radius(pos)
            cell_nrm /= r2

        num_photons_max = 10000000
        number_of_photons = np.zeros(num_cells, dtype="int64")
        energies = np.zeros(num_photons_max)

        start_e = 0
        end_e = 0

        spec = np.zeros(self.nbins)
        idxs = np.where(cut)[0]

        for ck in chunked(range(num_cells), 100):

            ibegin = ck[0]
            iend = ck[-1]+1
            nck = iend-ibegin

            cnm = cell_nrm[ibegin:iend]

            kTi = kT[ibegin:iend]

            if mode in ["photons", "spectrum"]:

                if self._density_dependence:
                    nHi = nH[ibegin:iend]
                    cspec, mspec, vspec = self.spectral_model.get_spectrum(kTi, nHi)
                else:
                    cspec, mspec, vspec = self.spectral_model.get_spectrum(kTi)

                tot_spec = cspec
                tot_spec += metalZ[ibegin:iend, np.newaxis] * mspec
                if self.num_var_elem > 0:
                    tot_spec += np.sum(elemZ[:, ibegin:iend, np.newaxis] * vspec, axis=0)

                if mode == "photons":

                    spec_sum = tot_spec.sum(axis=-1)
                    cell_norm = spec_sum * cnm

                    cell_n = np.atleast_1d(self.prng.poisson(lam=cell_norm))

                    number_of_photons[ibegin:iend] = cell_n
                    end_e += int(cell_n.sum())

                    norm_factor = 1.0 / spec_sum
                    p = norm_factor[:, np.newaxis]*tot_spec
                    cp = np.insert(np.cumsum(p, axis=-1), 0, 0.0, axis=1)
                    ei = start_e
                    for icell in range(nck):
                        cn = cell_n[icell]
                        if cn == 0:
                            continue
                        if self.method == "invert_cdf":
                            randvec = self.prng.uniform(size=cn)
                            randvec.sort()
                            cell_e = np.interp(randvec, cp[icell,:], self.bin_edges)
                        elif self.method == "accept_reject":
                            eidxs = self.prng.choice(self.nbins, size=cn, p=p[icell,:])
                            cell_e = self.emid[eidxs]
                        while ei+cn > num_photons_max:
                            num_photons_max *= 2
                        if num_photons_max > energies.size:
                            energies.resize(num_photons_max, refcheck=False)
                        energies[ei:ei+cn] = cell_e
                        ei += cn
                    start_e = end_e

                elif mode == "spectrum":

                    spec += np.sum(tot_spec * cnm[:, np.newaxis], axis=0)

            else:

                if self._density_dependence:
                    nHi = nH[ibegin:iend]
                    cflux, mflux, vflux = fluxf(kTi, nHi)
                else:
                    cflux, mflux, vflux = fluxf(kTi)

                tot_flux = cflux
                tot_flux += metalZ[ibegin:iend] * mflux
                if self.num_var_elem > 0:
                    tot_flux += np.sum(elemZ[:, ibegin:iend] * vflux, axis=0)

                ret[idxs[ibegin:iend]] = tot_flux*cnm

            self.pbar.update(nck)

        if mode == "photons":
            active_cells = number_of_photons > 0
            idxs = idxs[active_cells]
            ncells = idxs.size
            ee = energies[:end_e].copy()
            if self.binscale == "log":
                ee = 10**ee
            return ncells, number_of_photons[active_cells], idxs, ee
        elif mode == "spectrum":
            return spec
        else:
            return np.resize(ret, orig_shape)

    def cleanup_model(self):
        self.pbar.close()


class IGMSourceModel(ThermalSourceModel):
    _nei = False
    _density_dependence = True

    r"""
    A source model for a thermal plasma including photoionization and 
    resonant scattering from the CXB based on Khabibullin & Churazov 2019
    (https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4972K/) and Churazov 
    et al. 2001 (https://ui.adsabs.harvard.edu/abs/2001MNRAS.323...93C/).

    For temperatures higher than kT ~ 1.09 keV, a Cloudy-based CIE model
    is used to compute the spectrum. 

    Assumes the abundance tables from Feldman 1992.

    Table data and README files can be found at
    https://wwwmpa.mpa-garching.mpg.de/~ildar/igm/v3/.

    Parameters
    ----------
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The minimum energy for the spectral model.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
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
    nh_field : string or (ftype, fname) tuple
        The yt hydrogen nuclei density field (meaning all hydrogen, ionized or not) 
        to use for the model. Must have units of cm**-3. 
    resonant_scattering : boolean, optional
        Whether or not to include the effects of resonant scattering
        from CXB photons. Default: False
    cxb_factor : float, optional
        The fraction of the CXB photons that are resonant scattered to enhance
        the lines. Default: 0.5
    temperature_field : string or (ftype, fname) tuple, optional
        The yt temperature field to use for the thermal modeling. Must have
        units of Kelvin. If not specified, the default temperature field for
        the dataset will be used.
    emission_measure_field : string or (ftype, fname) tuple, optional
        The yt emission measure field to use for the thermal modeling. Must
        have units of cm^-3. If not specified, the default emission measure
        field for the dataset will be used or derived.
    h_fraction : float, string, or tuple of strings, optional
        The hydrogen mass fraction. If a float, assumes a constant mass 
        fraction of hydrogen throughout. If a string or tuple of strings, 
        is taken to be the name of the hydrogen fraction field. Defaults to
        the appropriate value for the Feldman abundance tables.
    kT_min : float, optional
        The default minimum temperature in keV to compute emission for.
        Default: 0.00431
    kT_max : float, optional
        The default maximum temperature in keV to compute emission for.
        Default: 64.0
    max_density : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The maximum density of the cells or particles to use when generating 
        photons. If a float, the units are assumed to be g/cm**3. 
        Default: None, meaning no maximum density.
    var_elem : dictionary, optional
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
    prng : integer or :class:`~numpy.random.RandomState` object 
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, 
        such as for a test. Default is to use the :mod:`numpy.random` module.
    """
    def __init__(self, emin, emax, nbins, Zmet, nh_field, binscale="linear",
                 resonant_scattering=False, cxb_factor=0.5, temperature_field=None,
                 emission_measure_field=None, h_fraction=None, kT_min=0.00431, 
                 kT_max=64.0, max_density=None, var_elem=None, method="invert_cdf",
                 prng=None):
        var_elem_keys = list(var_elem.keys()) if var_elem else None
        spectral_model = IGMSpectralModel(emin, emax, nbins, binscale=binscale,
                                          resonant_scattering=resonant_scattering,
                                          cxb_factor=cxb_factor, var_elem=var_elem_keys)
        nH_min = 10**spectral_model.Dvals[0]
        nH_max = 10**spectral_model.Dvals[-1]
        super().__init__(spectral_model, emin, emax, nbins, Zmet, binscale=binscale, kT_min=kT_min,
                         kT_max=kT_max, nH_min=nH_min, nH_max=nH_max, var_elem=var_elem,
                         max_density=max_density, method=method, abund_table="feld", prng=prng,
                         temperature_field=temperature_field, h_fraction=h_fraction,
                         emission_measure_field=emission_measure_field)
        self.nh_field = nh_field


class CIESourceModel(ThermalSourceModel):
    _nei = False
    _density_dependence = False
    r"""
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
        units of Kelvin. If not specified, the default temperature field for
        the dataset will be used.
    emission_measure_field : string or (ftype, fname) tuple, optional
        The yt emission measure field to use for the thermal modeling. Must
        have units of cm^-3. If not specified, the default emission measure
        field for the dataset will be used or derived.
    h_fraction : float, string, or tuple of strings, optional
        The hydrogen mass fraction. If a float, assumes a constant mass 
        fraction of hydrogen throughout. If a string or tuple of strings, 
        is taken to be the name of the hydrogen fraction field. Default is
        whatever value is appropriate for the chosen abundance tables.
    kT_min : float, optional
        The default minimum temperature in keV to compute emission for.
        Default: 0.025
    kT_max : float, optional
        The default maximum temperature in keV to compute emission for.
        Default: 64.0
    max_density : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The maximum density of the cells or particles to use when generating 
        photons. If a float, the units are assumed to be g/cm**3. 
        Default: None, meaning no maximum density.
    var_elem : dictionary, optional
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
        Whether or not the spectral lines should be thermally
        broadened. Only available for "apec" or "spex". Default: True
    model_root : string, optional
        The directory root where the model files are stored. If not provided,
        a default location known to pyXSIM is used.
    model_vers : string, optional
        The version identifier string for the model files, e.g.
        "2.0.2", if supported by the model. Currently only supported by
         "apec" and "spex". Default depends on the model being used.
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
    prng : integer or :class:`~numpy.random.RandomState` object 
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, 
        such as for a test. Default is to use the :mod:`numpy.random` module.

    Examples
    --------
    >>> source_model = CIESourceModel("apec", 0.1, 10.0, 10000,
    ...                               ("gas", "metallicity"))
    """
    def __init__(self, model, emin, emax, nbins, Zmet, binscale="linear", temperature_field=None,
                 emission_measure_field=None, h_fraction=None, kT_min=0.025,
                 kT_max=64.0, max_density=None, var_elem=None, method="invert_cdf",
                 thermal_broad=True, model_root=None, model_vers=None, nolines=False,
                 abund_table="angr", prng=None):
        var_elem_keys = list(var_elem.keys()) if var_elem else None
        if model in ["apec", "spex"]:
            spectral_model = TableCIEModel(model, emin, emax, nbins,
                                           binscale=binscale,
                                           var_elem=var_elem_keys,
                                           thermal_broad=thermal_broad,
                                           model_root=model_root,
                                           model_vers=model_vers,
                                           nolines=nolines, nei=self._nei,
                                           abund_table=abund_table)
        elif model == "mekal":
            spectral_model = MekalSpectralModel(emin, emax, nbins, binscale=binscale, 
                                                var_elem=var_elem_keys)
        elif model == "cloudy":
            if abund_table != "feld":
                mylog.warning("For the 'cloudy' model, the only available abundance table is "
                              "'feld', so using that one.")
                abund_table = "feld"
            spectral_model = CloudyCIESpectralModel(emin, emax, nbins, binscale=binscale,
                                                    var_elem=var_elem_keys)
        super().__init__(spectral_model, emin, emax, nbins, Zmet, binscale=binscale, kT_min=kT_min, 
                         kT_max=kT_max, var_elem=var_elem, max_density=max_density, method=method,
                         abund_table=abund_table, prng=prng, temperature_field=temperature_field,
                         emission_measure_field=emission_measure_field, h_fraction=h_fraction)
        self.var_elem_keys = self.spectral_model.var_elem_names
        self.var_ion_keys = self.spectral_model.var_ion_names


class NEISourceModel(CIESourceModel):
    _nei = True
    r"""
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
    var_elem : dictionary
        Abundances of elements. Each dictionary value, specified by the abundance 
        symbol, corresponds to the abundance of that symbol. If a float, it is 
        understood to be constant and in solar units. If a string or tuple of 
        strings, it is assumed to be a spatially varying field.
    binscale : string, optional
        The scale of the energy binning: "linear" or "log". 
        Default: "linear"
    temperature_field : string or (ftype, fname) tuple, optional
        The yt temperature field to use for the thermal modeling. Must have
        units of Kelvin. If not specified, the default temperature field for
        the dataset will be used.
    emission_measure_field : string or (ftype, fname) tuple, optional
        The yt emission measure field to use for the thermal modeling. Must
        have units of cm^-3. If not specified, the default emission measure
        field for the dataset will be used or derived.
    h_fraction : float, string, or tuple of strings, optional
        The hydrogen mass fraction. If a float, assumes a constant mass 
        fraction of hydrogen throughout. If a string or tuple of strings, 
        is taken to be the name of the hydrogen fraction field. Default is
        whatever value is appropriate for the chosen abundance tables.
    kT_min : float, optional
        The default minimum temperature in keV to compute emission for.
        Default: 0.025
    kT_max : float, optional
        The default maximum temperature in keV to compute emission for.
        Default: 64.0
    max_density : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The maximum density of the cells or particles to use when generating 
        photons. If a float, the units are assumed to be g/cm**3. 
        Default: None, meaning no maximum density.
    method : string, optional
        The method used to generate the photon energies from the spectrum:
        "invert_cdf": Invert the cumulative distribution function of the spectrum.
        "accept_reject": Acceptance-rejection method using the spectrum. 
        The first method should be sufficient for most cases.
    thermal_broad : boolean, optional
        Whether or not the spectral lines should be thermally
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
    >>> var_elem = {"H^1": ("flash", "h   "),
    >>>             "He^0": ("flash", "he  "),
    >>>             "He^1": ("flash", "he1 "),
    >>>             "He^2": ("flash", "he2 "),
    >>>             "O^0": ("flash", "o   "),
    >>>             "O^1": ("flash", "o1  "),
    >>>             "O^2": ("flash", "o2  "),
    >>>             "O^3": ("flash", "o3  "),
    >>>             "O^4": ("flash", "o4  "),
    >>>             "O^5": ("flash", "o5  "),
    >>>             "O^6": ("flash", "o6  "),
    >>>             "O^7": ("flash", "o7  "),
    >>>             "O^8": ("flash", "o8  ")
    >>>            }
    >>> source_model = ApecNEISourceModel(0.1, 10.0, 10000, var_elem)
    """
    def __init__(self, emin, emax, nbins, var_elem, binscale="linear", temperature_field=None,
                 emission_measure_field=None, h_fraction=None, kT_min=0.025,
                 kT_max=64.0, max_density=None, method="invert_cdf", thermal_broad=True,
                 model_root=None, model_vers=None, nolines=False, abund_table="angr", prng=None):
        super().__init__("apec", emin, emax, nbins, 0.0, binscale=binscale, 
                         temperature_field=temperature_field, 
                         emission_measure_field=emission_measure_field, h_fraction=h_fraction,
                         kT_min=kT_min, kT_max=kT_max, max_density=max_density, var_elem=var_elem,
                         method=method, thermal_broad=thermal_broad, model_root=model_root, 
                         model_vers=model_vers, nolines=nolines, abund_table=abund_table, 
                         prng=prng)
