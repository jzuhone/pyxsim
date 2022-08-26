from pyxsim.source_models.sources import SourceModel
from soxs.utils import parse_prng
from pyxsim.utils import parse_value
from yt.data_objects.static_output import Dataset

import numpy as np


class PowerLawSourceModel(SourceModel):
    r"""
    Initialize a source model from a power-law spectrum.

    Parameters
    ----------
    e0 : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The reference energy of the power law, in the rest frame of the source.
        If units are not given, they are assumed to be in keV.
    emin : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The minimum energy of the photons to be generated, in the rest frame of
        the source. If units are not given, they are assumed to be in keV.
    emax : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The maximum energy of the photons to be generated, in the rest frame of
        the source. If units are not given, they are assumed to be in keV.
    emission_field : string or (ftype, fname) tuple
        The field corresponding to the specific photon count rate per cell or
        particle, in the rest frame of the source, which serves as the
        normalization for the power law. Must be in counts/s/keV.
    index : float, string, or (ftype, fname) tuple
        The power-law index of the spectrum. Either a float for a single power law or
        the name of a field that corresponds to the power law.
    prng : integer or :class:`~numpy.random.RandomState` object 
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is to use the :mod:`numpy.random` module.

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
        self.prng = parse_prng(prng)
        self.spectral_norm = None
        self.redshift = None
        self.ftype = None

    def setup_model(self, data_source, redshift, spectral_norm):
        if isinstance(data_source, Dataset):
            ds = data_source
        else:
            ds = data_source.ds
        self.spectral_norm = spectral_norm
        self.redshift = redshift
        self.scale_factor = 1.0 / (1.0 + self.redshift)
        self.ftype = ds._get_field_info(self.emission_field).name[0]

    def make_spectrum(self, data_source, emin, emax, nbins, redshift=0.0,
                      dist=None, cosmology=None):
        ebins = np.linspace(emin, emax, nbins+1)
        spec = np.zeros(nbins)
        for chunk in data_source.chunks([], "io"):
            spec += self.process_data("spectrum", chunk, ebins=ebins)
        return self._make_spectrum(data_source.ds, ebins, spec,
                                   redshift, dist, cosmology)

    def make_fluxf(self, emin, emax, energy=False):
        return {"emin": emin, "emax": emax}

    def process_data(self, mode, chunk, fluxf=None, ebins=None):

        num_cells = len(chunk[self.emission_field])

        if isinstance(self.alpha, float):
            alpha = self.alpha*np.ones(num_cells)
        else:
            alpha = chunk[self.alpha].v

        if fluxf is None:
            ei = self.emin.v
            ef = self.emax.v
        else:
            ei = fluxf["emin"][0]
            ef = fluxf["emax"][1]

        if mode in ["photons", "photon_field"]:

            norm_fac = ef**(1.-alpha) - ei**(1.-alpha)
            norm_fac[alpha == 1] = np.log(ef / ei)
            norm_fac *= self.e0.v**alpha
            norm = norm_fac * chunk[self.emission_field].v
            norm[alpha != 1] /= (1.-alpha[alpha != 1])

            if mode == "photons":

                norm *= self.spectral_norm*self.scale_factor

                if self.observer == "internal":
                    pos = np.array([np.ravel(chunk[self.p_fields[i]].to_value("kpc"))
                                    for i in range(3)])
                    r2 = self.compute_radius(pos)
                    norm /= r2

                number_of_photons = self.prng.poisson(lam=norm)

                energies = np.zeros(number_of_photons.sum())

                start_e = 0
                end_e = 0
                for i in range(num_cells):
                    if number_of_photons[i] > 0:
                        end_e = start_e+number_of_photons[i]
                        u = self.prng.uniform(size=number_of_photons[i])
                        if alpha[i] == 1:
                            e = ei*(ef/ei)**u
                        else:
                            e = ei**(1.-alpha[i]) + u*norm_fac[i]
                            e **= 1./(1.-alpha[i])
                        energies[start_e:end_e] = e * self.scale_factor
                        start_e = end_e

                active_cells = number_of_photons > 0
                ncells = active_cells.sum()

                return ncells, number_of_photons[active_cells], active_cells, energies[:end_e].copy()

            elif mode == "photon_field":

                return norm

        elif mode == "energy_field":

            norm_fac = ef**(2.-alpha) - ei**(2.-alpha)
            norm_fac *= self.e0.v ** alpha / (2. - alpha)

            return norm_fac * chunk[self.emission_field].v

        elif mode == "spectrum":

            emid = 0.5*(ebins[1:]+ebins[:-1])
            return chunk[self.emission_field].v.sum()*(emid/self.e0.v)**(-alpha)
