"""
Classes for specific source models
"""
import numpy as np
from yt.funcs import get_pbar
from pyxsim.utils import mylog
from yt.units.yt_array import YTQuantity
from yt.utilities.physical_constants import mp, clight, kboltz
from pyxsim.utils import parse_value
from yt.utilities.exceptions import YTUnitConversionError

sqrt_two = np.sqrt(2.)

class SourceModel(object):

    def __init__(self, prng=None):
        self.spectral_norm = None
        self.redshift = None
        if prng is None:
            self.prng = np.random
        else:
            self.prng = prng

    def __call__(self, chunk):
        pass

    def setup_model(self, data_source, redshift, spectral_norm):
        self.spectral_norm = spectral_norm
        self.redshift = redshift

    def cleanup_model(self):
        self.spectral_norm = None
        self.redshift = None

particle_dens_fields = [("io", "density"),
                        ("PartType0", "Density")]
particle_temp_fields = [("io", "temperature"),
                        ("PartType0", "Temperature")]

class ThermalSourceModel(SourceModel):
    r"""
    Initialize a source model from a thermal spectrum.

    Parameters
    ----------
    spectral_model : :class:`~pyxsim.spectral_models.SpectralModel`
        A thermal spectral model instance, either of :class:`~pyxsim.spectral_models.XSpecThermalModel` or :class:`~pyxsim.spectral_models.TableApecModel`.
    temperature_field : string or (ftype, fname) tuple, optional
        The yt temperature field to use for the thermal modeling. Must have units
        of Kelvin. If not specified, the default temperature field for the dataset
        will be used.
    emission_measure_field : string or (ftype, fname) tuple, optional
        The yt emission measure field to use for the thermal modeling. Must have units
        of cm^-3. If not specified, the default emission measure field for the dataset
        will be used or derived.
    kT_min : float, optional
        The default minimum temperature in keV to compute emission for. Default: 0.008
    kT_max : float, optional
        The default maximum temperature in keV to compute emission for. Default: 64.0
    n_kT : integer, optional
        The number of temperature bins to use when computing emission. Default: 10000
    kT_scale : string, optional
        The scaling of the bins to use when computing emission, "linear" or "log". 
        Default: "linear"
    Zmet : float or string, optional
        The metallicity. If a float, assumes a constant metallicity throughout.
        If a string, is taken to be the name of the metallicity field.
    method : string, optional
        The method used to generate the photon energies from the spectrum:
        "invert_cdf": Invert the cumulative distribution function of the spectrum.
        "accept_reject": Acceptance-rejection method using the spectrum. 
        The first method should be sufficient for most cases. 
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.

    Examples
    --------
    >>> mekal_model = XSpecThermalModel("mekal", 0.05, 50.0, 1000)
    >>> source_model = ThermalSourceModel(mekal_model, Zmet="metallicity")
    """
    def __init__(self, spectral_model, temperature_field=None,
                 emission_measure_field=None, kT_min=0.008,
                 kT_max=64.0, n_kT=10000, kT_scale="linear", 
                 Zmet=0.3, method="invert_cdf", prng=None):
        self.temperature_field = temperature_field
        self.Zmet = Zmet
        self.spectral_model = spectral_model
        self.method = method
        if prng is None:
            self.prng = np.random
        else:
            self.prng = prng
        self.kT_min = kT_min
        self.kT_max = kT_max
        self.kT_scale = kT_scale
        self.n_kT = n_kT
        self.spectral_norm = None
        self.redshift = None
        self.pbar = None
        self.kT_bins = None
        self.dkT = None
        self.emission_measure_field = emission_measure_field

    def setup_model(self, data_source, redshift, spectral_norm):
        self.redshift = redshift
        if self.emission_measure_field is None:
            found_dfield = [fd for fd in particle_dens_fields if fd in data_source.ds.field_list]
            if len(found_dfield) > 0:
                ptype = found_dfield[0][0]
                def _emission_measure(field, data):
                    nenh = data[found_dfield[0]]*data['particle_mass']
                    nenh /= mp*mp
                    nenh.convert_to_units("cm**-3")
                    if data.has_field_parameter("X_H"):
                        X_H = data.get_field_parameter("X_H")
                    else:
                        X_H = 0.76
                    if ('PartType0', 'ElectronAbundance') in data_source.ds.field_list:
                        nenh *= X_H * data['PartType0','ElectronAbundance']
                        nenh *= X_H * (1.-data['PartType0','NeutralHydrogenAbundance'])
                    else:
                        nenh *= 0.5*(1.+X_H)*X_H
                    return nenh
                data_source.ds.add_field((ptype, 'emission_measure'),
                                         function=_emission_measure,
                                         particle_type=True,
                                         units="cm**-3")
                self.emission_measure_field = (ptype, 'emission_measure')
            else:
                self.emission_measure_field = ('gas', 'emission_measure')
        mylog.info("Using emission measure field '(%s, %s)'." % self.emission_measure_field)
        if self.temperature_field is None:
            found_tfield = [fd for fd in particle_temp_fields if fd in data_source.ds.derived_field_list]
            if len(found_tfield) > 0:
                self.temperature_field = found_tfield[0]
            else:
                self.temperature_field = ('gas', 'temperature')
        mylog.info("Using temperature field '(%s, %s)'." % self.temperature_field)
        self.spectral_model.prepare_spectrum(redshift)
        self.spectral_norm = spectral_norm
        if self.kT_scale == "linear":
            self.kT_bins = np.linspace(self.kT_min, self.kT_max, num=self.n_kT+1)
        elif self.kT_scale == "log":
            self.kT_bins = np.logspace(np.log10(self.kT_min), np.log10(self.kT_max), 
                                       num=self.n_kT+1)
        self.dkT = np.diff(self.kT_bins)
        kT = (kboltz*data_source[self.temperature_field]).in_units("keV").v
        num_cells = np.logical_and(kT > self.kT_min, kT < self.kT_max).sum()
        self.source_type = data_source.ds._get_field_info(self.emission_measure_field).name[0]
        self.pbar = get_pbar("Generating photons ", num_cells)

    def __call__(self, chunk):

        num_photons_max = 10000000
        emid = self.spectral_model.emid
        ebins = self.spectral_model.ebins
        nchan = len(emid)

        kT = (kboltz*chunk[self.temperature_field]).in_units("keV").v
        if len(kT) == 0:
            return
        EM = chunk[self.emission_measure_field].v

        idxs = np.argsort(kT)

        kT_sorted = kT[idxs]
        idx_min = np.searchsorted(kT_sorted, self.kT_min)
        idx_max = np.searchsorted(kT_sorted, self.kT_max)
        idxs = idxs[idx_min:idx_max]
        num_cells = len(idxs)
        if num_cells == 0:
            return

        kT_idxs = np.digitize(kT[idxs], self.kT_bins)-1
        bcounts = np.bincount(kT_idxs).astype("int")
        bcounts = bcounts[bcounts > 0]
        n = int(0)
        bcell = []
        ecell = []
        for bcount in bcounts:
            bcell.append(n)
            ecell.append(n+bcount)
            n += bcount
        kT_idxs = np.unique(kT_idxs)

        cell_em = EM[idxs]*self.spectral_norm

        if isinstance(self.Zmet, float):
            metalZ = self.Zmet*np.ones(num_cells)
        else:
            metalZ = chunk[self.Zmet].v[idxs]

        number_of_photons = np.zeros(num_cells, dtype="uint64")
        energies = np.zeros(num_photons_max)

        start_e = 0
        end_e = 0

        for ibegin, iend, ikT in zip(bcell, ecell, kT_idxs):

            kT = self.kT_bins[ikT] + 0.5*self.dkT[ikT]

            n_current = iend-ibegin

            cem = cell_em[ibegin:iend]

            cspec, mspec = self.spectral_model.get_spectrum(kT)

            tot_ph_c = cspec.d.sum()
            tot_ph_m = mspec.d.sum()

            u = self.prng.uniform(size=n_current)

            cell_norm_c = tot_ph_c*cem
            cell_norm_m = tot_ph_m*metalZ[ibegin:iend]*cem
            cell_norm = np.modf(cell_norm_c + cell_norm_m)
            cell_n = np.uint64(cell_norm[1]) + np.uint64(cell_norm[0] >= u)

            number_of_photons[ibegin:iend] = cell_n

            end_e += int(cell_n.sum())

            if self.method == "invert_cdf":
                cumspec_c = np.cumsum(cspec.d)
                cumspec_m = np.cumsum(mspec.d)
                cumspec_c = np.insert(cumspec_c, 0, 0.0)
                cumspec_m = np.insert(cumspec_m, 0, 0.0)

            ei = start_e
            for cn, Z in zip(number_of_photons[ibegin:iend], metalZ[ibegin:iend]):
                self.pbar.update()
                if cn == 0:
                    continue
                # The rather verbose form of the few next statements is a
                # result of code optimization and shouldn't be changed
                # without checking for perfomance degradation. See
                # https://bitbucket.org/yt_analysis/yt/pull-requests/1766
                # for details.
                if self.method == "invert_cdf":
                    cumspec = cumspec_c
                    cumspec += Z * cumspec_m
                    norm_factor = 1.0 / cumspec[-1]
                    cumspec *= norm_factor
                    randvec = self.prng.uniform(size=cn)
                    randvec.sort()
                    cell_e = np.interp(randvec, cumspec, ebins)
                elif self.method == "accept_reject":
                    tot_spec = cspec.d
                    tot_spec += Z * mspec.d
                    norm_factor = 1.0 / tot_spec.sum()
                    tot_spec *= norm_factor
                    eidxs = self.prng.choice(nchan, size=cn, p=tot_spec)
                    cell_e = emid[eidxs]
                while ei+cn > num_photons_max:
                    num_photons_max *= 2
                if num_photons_max > energies.size:
                    energies.resize(num_photons_max, refcheck=False)
                energies[ei:ei+cn] = cell_e
                ei += cn

            start_e = end_e

        active_cells = number_of_photons > 0
        idxs = idxs[active_cells]

        return number_of_photons[active_cells], idxs, energies[:end_e].copy()

    def cleanup_model(self):
        self.pbar.finish()
        self.redshift = None
        self.spectral_model.cleanup_spectrum()
        self.pbar = None
        self.spectral_norm = None
        self.kT_bins = None
        self.dkT = None

class PowerLawSourceModel(SourceModel):
    r"""
    Initialize a source model from a power-law spectrum.

    Parameters
    ----------
    e0 : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
        The reference energy of the power law. If units are not given,
        they are assumed to be in keV.
    emin : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
        The minimum energy of the photons to be generated. If units
        are not given, they are assumed to be in keV.
    emax : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
        The maximum energy of the photons to be generated. If units
        are not given, they are assumed to be in keV.
    emission_field : string or (ftype, fname) tuple
        The field which serves as the normalization for the power law. Must be in units
        of counts/s/keV.
    index : float, string, or (ftype, fname) tuple
        The power-law index of the spectrum. Either a float for a single power law or
        the name of a field that corresponds to the power law.
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.

    Examples
    --------
    >>> e0 = (1.0, "keV")
    >>> emin = (0.01, "keV")
    >>> emax = (100., "keV")
    >>> plaw_model = PowerLawSourceModel(e0, emin, emax, ("gas", "norm"), ("gas", "index"))
    """
    def __init__(self, e0, emin, emax, emission_field, alpha, prng=None):
        self.e0 = parse_value(e0, "keV")
        self.emin = parse_value(emin, "keV")
        self.emax = parse_value(emax, "keV")
        self.emission_field = emission_field
        self.alpha = alpha
        if prng is None:
            self.prng = np.random
        else:
            self.prng = prng
        self.spectral_norm = None
        self.redshift = None

    def setup_model(self, data_source, redshift, spectral_norm):
        self.spectral_norm = spectral_norm
        self.redshift = redshift
        self.source_type = data_source.ds._get_field_info(self.emission_field).name[0]

    def __call__(self, chunk):

        num_cells = len(chunk[self.emission_field])

        if isinstance(self.alpha, float):
            alpha = self.alpha*np.ones(num_cells)
        else:
            alpha = chunk[self.alpha].v

        norm_fac = (self.emax.v**(1.-alpha)-self.emin.v**(1.-alpha))
        norm_fac[alpha == 1] = np.log(self.emax.v/self.emin.v)
        norm = norm_fac*chunk[self.emission_field].v*self.e0.v**alpha
        norm[alpha != 1] /= (1.-alpha[alpha != 1])
        norm *= self.spectral_norm
        norm = np.modf(norm)

        u = self.prng.uniform(size=num_cells)
        number_of_photons = np.uint64(norm[1]) + np.uint64(norm[0] >= u)

        energies = np.zeros(number_of_photons.sum())

        start_e = 0
        end_e = 0
        for i in range(num_cells):
            if number_of_photons[i] > 0:
                end_e = start_e+number_of_photons[i]
                u = self.prng.uniform(size=number_of_photons[i])
                if alpha[i] == 1:
                    e = self.emin.v*(self.emax.v/self.emin.v)**u
                else:
                    e = self.emin.v**(1.-alpha[i]) + u*norm_fac[i]
                    e **= 1./(1.-alpha[i])
                energies[start_e:end_e] = e / (1.+self.redshift)
                start_e = end_e

        active_cells = number_of_photons > 0

        return number_of_photons[active_cells], active_cells, energies[:end_e].copy()

    def cleanup_model(self):
        self.redshift = None
        self.spectral_norm = None

class LineSourceModel(SourceModel):
    r"""
    Initialize a source model from a single line.

    Parameters
    ----------
    e0 : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
        The location of the emission line in energy in the rest frame of the
        source. If units are not given, they are assumed to be in keV.
    emission_field : string or (ftype, fname) tuple
        The field which serves as the normalization for the line. Must be in
        counts/s.
    sigma : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or field name, optional
        The standard intrinsic deviation of the emission line (not from Doppler
        broadening, which is handled in the projection step). Units of
        velocity or energy are accepted. If units are not given, they
        are assumed to be in keV. If set to a field name, the line broadening
        is assumed to be based on this field (in units of velocity or energy).
        If set to None (the default), it is assumed that the line is unbroadened.
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.

    Examples
    --------
    >>> location = (3.5, "keV")
    >>> sigma = (1000., "km/s")
    >>> line_model = LineEmissionSourceModel(location, "dark_matter_density_squared", sigma=sigma)
    """
    def __init__(self, e0, emission_field, sigma=None, prng=None):
        self.e0 = parse_value(e0, "keV")
        if isinstance(sigma, (float, YTQuantity)) or (isinstance(sigma, tuple) and isinstance(sigma[0], float)):
            # The broadening is constant
            try:
                self.sigma = parse_value(sigma, "keV")
            except YTUnitConversionError:
                try:
                    self.sigma = parse_value(sigma, "km/s")
                    self.sigma *= self.e0/clight
                    self.sigma.convert_to_units("keV")
                except YTUnitConversionError:
                    raise RuntimeError("Units for sigma must either be in dimensions of "
                                       "energy or velocity! sigma = %s" % sigma)
        else:
            # Either no broadening or a field name
            self.sigma = sigma
        self.emission_field = emission_field
        if prng is None:
            self.prng = np.random
        else:
            self.prng = prng
        self.spectral_norm = None
        self.redshift = None

    def setup_model(self, data_source, redshift, spectral_norm):
        self.spectral_norm = spectral_norm
        self.redshift = redshift
        self.source_type = data_source.ds._get_field_info(self.emission_field).name[0]

    def __call__(self, chunk):
        num_cells = len(chunk[self.emission_field])
        F = chunk[self.emission_field]*self.spectral_norm
        norm = np.modf(F.in_cgs().v)
        u = self.prng.uniform(size=num_cells)
        number_of_photons = np.uint64(norm[1]) + np.uint64(norm[0] >= u)

        energies = self.e0*np.ones(number_of_photons.sum())

        if isinstance(self.sigma, YTQuantity):
            dE = self.prng.normal(loc=0.0, scale=float(self.sigma),
                                  size=number_of_photons.sum())*self.e0.uq
            energies += dE
        elif self.sigma is not None:
            sigma = (chunk[self.sigma]*self.e0/clight).in_units("keV")
            start_e = 0
            for i in range(num_cells):
                if number_of_photons[i] > 0:
                    end_e = start_e+number_of_photons[i]
                    dE = self.prng.normal(loc=0.0, scale=float(sigma[i]),
                                          size=number_of_photons[i])*self.e0.uq
                    energies[start_e:end_e] += dE
                    start_e = end_e

        energies = energies / (1.+self.redshift)

        active_cells = number_of_photons > 0

        return number_of_photons[active_cells], active_cells, energies

    def cleanup_model(self):
        self.redshift = None
        self.spectral_norm = None
