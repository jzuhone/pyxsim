"""
Classes for specific photon models

The algorithms used here are based off of the method used by the
PHOX code (http://www.mpa-garching.mpg.de/~kdolag/Phox/),
developed by Veronica Biffi and Klaus Dolag. References for
PHOX may be found at:

Biffi, V., Dolag, K., Bohringer, H., & Lemson, G. 2012, MNRAS, 420, 3545
http://adsabs.harvard.edu/abs/2012MNRAS.420.3545B

Biffi, V., Dolag, K., Bohringer, H. 2013, MNRAS, 428, 1395
http://adsabs.harvard.edu/abs/2013MNRAS.428.1395B

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from six import string_types
import numpy as np
from yt.units.yt_array import YTQuantity
from yt.utilities.physical_constants import mp, clight, kboltz
from yt.analysis_modules.photon_simulator.photon_simulator import \
    parse_value
from yt.utilities.exceptions import YTUnitConversionError

n_kT = 10000
kT_min = 8.08e-2
kT_max = 50.
sqrt_two = np.sqrt(2.)

class SourceModel(object):

    def __init__(self):
        self.spectral_norm = None
        self.pbar = None

    def __call__(self, chunk):
        pass

    def setup_model(self, data_source, redshift, spectral_norm, pbar):
        self.spectral_norm = spectral_norm
        self.pbar = pbar

    def cleanup_model(self):
        self.spectral_norm = None
        self.pbar = None

class ThermalSourceModel(SourceModel):
    r"""
    Initialize a ThermalSourceModel from a thermal spectrum.

    Parameters
    ----------
    spectral_model : `SpectralModel`
        A thermal spectral model instance, either of `XSpecThermalModel`
        or `TableApecModel`.
    Zmet : float or string, optional
        The metallicity. If a float, assumes a constant metallicity throughout.
        If a string, is taken to be the name of the metallicity field.
    photons_per_chunk : integer
        The maximum number of photons that are allocated per chunk. Increase or decrease
        as needed.
    method : string, optional
        The method used to generate the photon energies from the spectrum:
        "invert_cdf": Invert the cumulative distribution function of the spectrum.
        "accept_reject": Acceptance-rejection method using the spectrum. 
        The first method should be sufficient for most cases. 
    prng : NumPy `RandomState` object or numpy.random
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the numpy.random module.
    velocity_fields : list of field names
        The fields to use from the data object as the velocity fields for Doppler
        broadening. Default is ["velocity_x","velocity_y","velocity_z"].

    Examples
    --------
    >>> mekal_model = XSpecThermalModel("mekal", 0.05, 50.0, 1000)
    >>> source_model = ThermalSourceModel(mekal_model, X_H=0.76,
    ...                                   Zmet="metallicity")
    """
    def __init__(self, spectral_model, temperature_field="temperature",
                 emission_measure_field="emission_measure",
                 Zmet=0.3, photons_per_chunk=10000000,
                 method="invert_cdf", prng=np.random):
        self.temperature_field = temperature_field
        self.Zmet = Zmet
        self.spectral_model = spectral_model
        self.photons_per_chunk = photons_per_chunk
        self.method = method
        self.prng = prng
        self.spectral_norm = None
        self.pbar = None
        self.kT_bins = None
        self.dkT = None
        self.emission_measure_field = emission_measure_field

    def setup_model(self, data_source, redshift, spectral_norm, pbar):
        my_kT_min, my_kT_max = data_source.quantities.extrema("kT")
        self.spectral_model.prepare_spectrum(redshift)
        self.spectral_norm = spectral_norm
        self.pbar = pbar
        self.kT_bins = np.linspace(kT_min, max(my_kT_max.v, kT_max), num=n_kT+1)
        self.dkT = self.kT_bins[1]-self.kT_bins[0]

    def __call__(self, chunk):

        emid = self.spectral_model.emid
        ebins = self.spectral_model.ebins
        nchan = len(emid)

        kT = (kboltz*chunk[self.temperature_field]).in_units("keV").v
        num_cells = len(kT)
        if num_cells == 0:
            return
        EM = chunk[self.emission_measure_field].v

        if isinstance(self.Zmet, string_types):
            metalZ = chunk[self.Zmet].v
        else:
            metalZ = self.Zmet*np.ones(num_cells)

        idxs = np.argsort(kT)

        kT_idxs = np.digitize(kT[idxs], self.kT_bins)
        kT_idxs = np.minimum(np.maximum(1, kT_idxs), n_kT) - 1
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

        number_of_photons = np.zeros(num_cells, dtype="uint64")
        energies = np.zeros(self.photons_per_chunk)

        start_e = 0
        end_e = 0

        for ibegin, iend, ikT in zip(bcell, ecell, kT_idxs):

            kT = self.kT_bins[ikT] + 0.5*self.dkT

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

            if end_e > self.photons_per_chunk:
                raise RuntimeError("Number of photons generated for this chunk "+
                                   "exceeds photons_per_chunk (%d)! " % self.photons_per_chunk +
                                   "Increase photons_per_chunk!")

            if self.method == "invert_cdf":
                cumspec_c = np.cumsum(cspec.d)
                cumspec_m = np.cumsum(mspec.d)
                cumspec_c = np.insert(cumspec_c, 0, 0.0)
                cumspec_m = np.insert(cumspec_m, 0, 0.0)

            ei = start_e
            for cn, Z in zip(number_of_photons[ibegin:iend], metalZ[ibegin:iend]):
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
                energies[ei:ei+cn] = cell_e
                self.pbar.update()
                ei += cn

            start_e = end_e

        active_cells = number_of_photons > 0

        return number_of_photons[active_cells], active_cells, energies[:end_e].copy()

    def cleanup_model(self):
        self.spectral_model.cleanup_spectrum()
        self.pbar = None
        self.spectral_norm = None
        self.kT_bins = None
        self.dkT = None

class PowerLawSourceModel(SourceModel):
    r"""
    Initialize a PowerLawSourceModel from a power-law spectrum.

    Parameters
    ----------
    e0 : float, (value, unit) tuple, or YTQuantity
        The reference energy of the power law. If units are not given,
        they are assumed to be in keV.
    emin : float, (value, unit) tuple, or YTQuantity
        The minimum energy of the photons to be generated. If units
        are not given, they are assumed to be in keV.
    emax : float, (value, unit) tuple, or YTQuantity
        The maximum energy of the photons to be generated. If units
        are not given, they are assumed to be in keV.
    norm_field : string or (ftype, fname) tuple
        The field which serves as the normalization for the power law. Must be in units
        of counts/s/cm**3/keV.
    index : float, string, or (ftype, fname) tuple
        The power-law index of the spectrum. Either a float for a single power law or
        the name of a field that corresponds to the power law.
    prng : NumPy `RandomState` object or numpy.random
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is the numpy.random module.

    Examples
    --------
    >>> e0 = (1.0, "keV")
    >>> emin = (0.01, "keV")
    >>> emax = (100., "keV")
    >>> plaw_model = PowerLawSourceModel(e0, emin, emax, ("gas", "norm"), ("gas", "index"))
    """
    def __init__(self, e0, emin, emax, norm_field, index, prng=np.random):
        self.e0 = parse_value(e0, "keV")
        self.emin = parse_value(emin, "keV")
        self.emax = parse_value(emax, "keV")
        self.norm_field = norm_field
        self.index = index
        self.prng = prng
        self.spectral_norm = None
        self.pbar = None
        self.redshift = None

    def setup_model(self, data_source, redshift, spectral_norm, pbar):
        self.spectral_norm = spectral_norm
        self.redshift = redshift
        self.pbar = pbar

    def __call__(self, chunk):

        if isinstance(self.index, float):
            index = self.index
        else:
            index = chunk[self.index].v
        norm_fac = (self.emax**(1.-index)-self.emin**(1.-index))*self.spectral_norm
        num_cells = len(norm_fac)

        norm = norm_fac*chunk[self.norm_field]*self.e0**index/(1.-index)
        norm = np.modf(norm.in_cgs().v)
        u = self.prng.uniform(size=num_cells)
        number_of_photons = np.uint64(norm[1]) + np.uint64(norm[0] >= u)

        energies = np.zeros(number_of_photons.sum())

        start_e = 0
        for i in range(num_cells):
            if number_of_photons[i] > 0:
                end_e = start_e+number_of_photons[i]
                u = self.prng.uniform(size=number_of_photons[i])
                e = self.emin**(1.-index[i]) + u*norm[i]
                e **= 1./(1.-index[i])
                energies[start_e:end_e] = e / (1.+self.redshift)
                start_e = end_e
            self.pbar.update()

        active_cells = number_of_photons > 0

        return number_of_photons[active_cells], active_cells, energies[:end_e].copy()

    def cleanup_model(self):
        self.redshift = None
        self.pbar = None
        self.spectral_norm = None

class LineEmissionSourceModel(SourceModel):
    r"""
    Initialize a LineEmissionSourceModel from a single line.

    Parameters
    ----------
    location : float, (value, unit) tuple, or YTQuantity
        The location of the emission line in energy in the rest frame of the
        object. If units are not given, they are assumed to be in keV.
    amplitude_field : string or (ftype, fname) tuple
        The field which serves as the normalization for the linej. Must be in
        counts/s/cm**3.
    sigma : float, (value, unit) tuple, YTQuantity, or field name, optional
        The standard intrinsic deviation of the emission line (not from Doppler
        broadening, which is handled in the projection step). Units of
        velocity or energy are accepted. If units are not given, they
        are assumed to be in keV. If set to a field name, the line broadening
        is assumed to be based on this field (in units of velocity or energy).
        If set to None (the default), it is assumed that the line is unbroadened.
    prng : NumPy `RandomState` object or numpy.random
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is the numpy.random module.

    Examples
    --------
    >>> location = (3.5, "keV")
    >>> sigma = (1000., "km/s")
    >>> line_model = LineEmissionSourceModel(location, "dark_matter_density_squared", sigma=sigma)
    """
    def __init__(self, location, amplitude_field, sigma=None, prng=np.random):
        self.location = parse_value(location, "keV")
        if isinstance(sigma, (float, YTQuantity)) or (isinstance(sigma, tuple) and isinstance(sigma[0], float)):
            # The broadening is constant
            try:
                self.sigma = parse_value(sigma, "keV")
            except YTUnitConversionError:
                try:
                    self.sigma = parse_value(sigma, "km/s")
                    self.sigma *= self.location/clight
                    self.sigma.convert_to_units("keV")
                except YTUnitConversionError:
                    raise RuntimeError("Units for sigma must either be in dimensions of "
                                       "energy or velocity! sigma = %s" % sigma)
        else:
            # Either no broadening or a field name
            self.sigma = sigma
        self.amplitude_field = amplitude_field
        self.prng = prng
        self.spectral_norm = None
        self.pbar = None
        self.redshift = None

    def setup_model(self, data_source, redshift, spectral_norm, pbar):
        self.spectral_norm = spectral_norm
        self.redshift = redshift
        self.pbar = pbar

    def __call__(self, chunk):
        num_cells = len(chunk["x"])
        F = chunk[self.amplitude_field]*chunk["cell_volume"]*self.spectral_norm
        norm = np.modf(F.in_cgs().v)
        u = self.prng.uniform(size=num_cells)
        number_of_photons = np.uint64(norm[1]) + np.uint64(norm[0] >= u)

        energies = self.location*np.ones(number_of_photons.sum())

        if isinstance(self.sigma, YTQuantity):
            dE = self.prng.normal(loc=0.0, scale=self.sigma.v,
                                  size=number_of_photons.sum())*self.location.uq
            energies += dE
        elif self.sigma is not None:
            start_e = 0
            for i in range(num_cells):
                if number_of_photons[i] > 0:
                    end_e = start_e+number_of_photons[i]
                    dE = self.prng.normal(loc=0.0, scale=chunk[self.sigma][i].v,
                                          size=number_of_photons[i])*self.location.uq
                    energies[start_e:end_e] += dE
                    start_e = end_e
            self.pbar.update()

        energies = energies / (1.+self.redshift)

        active_cells = number_of_photons > 0

        return number_of_photons[active_cells], active_cells, energies

    def cleanup_model(self):
        self.redshift = None
        self.pbar = None
        self.spectral_norm = None