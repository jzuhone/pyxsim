from numbers import Number

import numpy as np
from soxs.utils import parse_prng
from yt.data_objects.static_output import Dataset

from pyxsim.lib.spectra import power_law_spectrum
from pyxsim.source_models.sources import SourceModel
from pyxsim.utils import mylog, parse_value


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
        self.ftype = None

    def setup_model(self, mode, data_source, redshift):
        if isinstance(data_source, Dataset):
            ds = data_source
        else:
            ds = data_source.ds
        self.scale_factor = 1.0 / (1.0 + redshift)
        self.emission_field = ds._get_field_info(self.emission_field).name
        if not isinstance(self.alpha, Number):
            self.alpha = ds._get_field_info(self.alpha).name
            if self.emission_field[0] != self.alpha[0]:
                mylog.warning(
                    "The 'emission_field' %s and the 'alpha' field %s do not have the same field type!",
                    self.emission_field,
                    self.alpha,
                )
        self.ftype = self.emission_field[0]
        if mode == "spectrum":
            self.setup_pbar(data_source, self.emission_field)

    def __repr__(self):
        rets = [
            "PowerLawSourceModel(\n",
            f"    e0={self.e0}\n",
            f"    emin={self.emin}\n",
            f"    emax={self.emax}\n",
            f"    emission_field={self.emission_field}\n",
            f"    alpha={self.alpha}\n",
            ")",
        ]
        return "".join(rets)

    def cleanup_model(self, mode):
        if mode == "spectrum":
            self.pbar.close()

    def make_spectrum(
        self, data_source, emin, emax, nbins, redshift=0.0, dist=None, cosmology=None
    ):
        """
        Make a count rate spectrum in the source frame from a yt data container,
        or a spectrum in the observer frame.

        Parameters
        ----------
        data_source : :class:`~yt.data_objects.data_containers.YTSelectionContainer`
            The data source from which the photons will be generated.
        emin : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The minimum energy in the band. If a float, it is assumed to be
            in keV.
        emax : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The minimum energy in the band. If a float, it is assumed to be
            in keV.
        nbins : integer
            The number of bins in the spectrum.
        redshift : float, optional
            If greater than 0, we assume that the spectrum should be created in
            the observer frame at a distance given by the cosmology. Default: 0.0
        dist : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity, optional
            The distance to a nearby source, if redshift = 0.0. If a float, it
            is assumed to be in units of kpc.
        cosmology : :class:`~yt.utilities.cosmology.Cosmology`, optional
            Cosmological information. If not supplied, we try to get the
            cosmology from the dataset. Otherwise, LCDM with the default yt
            parameters is assumed.

        Returns
        -------
        :class:`~soxs.spectra.CountRateSpectrum` or :class:`~soxs.spectra.Spectrum`,
        depending on how the method is invoked.
        """
        ebins = np.linspace(emin, emax, nbins + 1)
        spec = np.zeros(nbins)
        spectral_norm = 1.0
        self.setup_model("spectrum", data_source, redshift)
        for chunk in data_source.chunks([], "io"):
            spec += self.process_data("spectrum", chunk, spectral_norm, ebins=ebins)
        self.cleanup_model("spectrum")
        return self._make_spectrum(
            data_source.ds, ebins, spec, redshift, dist, cosmology
        )

    def make_fluxf(self, emin, emax, energy=False):
        return {"emin": emin, "emax": emax}

    def process_data(self, mode, chunk, spectral_norm, fluxf=None, ebins=None):

        num_cells = len(chunk[self.emission_field])

        if isinstance(self.alpha, float):
            alpha = self.alpha * np.ones(num_cells)
        else:
            alpha = chunk[self.alpha].d

        if fluxf is None:
            ei = self.emin.v
            ef = self.emax.v
        else:
            ei = fluxf["emin"].v
            ef = fluxf["emax"].v

        if mode in ["photons", "photon_field"]:

            norm_fac = ef ** (1.0 - alpha) - ei ** (1.0 - alpha)
            norm_fac[alpha == 1] = np.log(ef / ei)
            norm_fac *= self.e0.v**alpha
            norm = norm_fac * chunk[self.emission_field].d
            if np.any(alpha != 1):
                norm[alpha != 1] /= 1.0 - alpha[alpha != 1]

            if mode == "photons":

                norm *= spectral_norm * self.scale_factor

                if self.observer == "internal":
                    pos = np.array(
                        [
                            np.ravel(chunk[self.p_fields[i]].to_value("kpc"))
                            for i in range(3)
                        ]
                    )
                    r2 = self.compute_radius(pos)
                    norm /= r2

                number_of_photons = self.prng.poisson(lam=norm)

                energies = np.zeros(number_of_photons.sum())

                start_e = 0
                end_e = 0
                for i in range(num_cells):
                    if number_of_photons[i] > 0:
                        end_e = start_e + number_of_photons[i]
                        u = self.prng.uniform(size=number_of_photons[i])
                        if alpha[i] == 1:
                            e = ei * (ef / ei) ** u
                        else:
                            e = ei ** (1.0 - alpha[i]) + u * norm_fac[i]
                            e **= 1.0 / (1.0 - alpha[i])
                        energies[start_e:end_e] = e * self.scale_factor
                        start_e = end_e

                active_cells = number_of_photons > 0
                ncells = active_cells.sum()

                return (
                    ncells,
                    number_of_photons[active_cells],
                    active_cells,
                    energies[:end_e].copy(),
                )

            elif mode == "photon_field":

                return norm

        elif mode == "energy_field":

            norm_fac = ef ** (2.0 - alpha) - ei ** (2.0 - alpha)
            norm_fac *= self.e0.v**alpha / (2.0 - alpha)

            return norm_fac * chunk[self.emission_field].d

        elif mode == "spectrum":

            inv_sf = 1.0 / self.scale_factor
            emid = 0.5 * (ebins[1:] + ebins[:-1]) * inv_sf / self.e0.v

            spec = power_law_spectrum(
                num_cells, emid, alpha, chunk[self.emission_field].d, self.pbar
            )

            return spec
