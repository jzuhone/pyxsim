from numbers import Number

import numpy as np
from soxs.utils import parse_prng
from unyt.exceptions import UnitConversionError
from yt.data_objects.static_output import Dataset

from pyxsim.lib.spectra import power_law_spectrum
from pyxsim.source_models.sources import SourceModel
from pyxsim.utils import check_num_cells, mylog, parse_value, sanitize_normal


class PowerLawSourceModel(SourceModel):
    """
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
    luminosity_field : string or (ftype, fname) tuple
        The field corresponding to the luminosity within the emin-emax band per
        cell or particle, in the rest frame of the source, which serves as the
        normalization for the power law. Must be in units with dimensions of power,
        such as erg/s, W, or keV/s.
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

    def __init__(self, e0, emin, emax, luminosity_field, alpha, prng=None):
        self.e0 = parse_value(e0, "keV")
        self.emin = parse_value(emin, "keV")
        self.emax = parse_value(emax, "keV")
        self.luminosity_field = luminosity_field
        self.alpha = alpha
        self.prng = parse_prng(prng)
        self.ftype = None

    def setup_model(self, mode, data_source, redshift):
        if isinstance(data_source, Dataset):
            ds = data_source
        else:
            ds = data_source.ds
        self.scale_factor = 1.0 / (1.0 + redshift)
        self.luminosity_field = ds._get_field_info(self.luminosity_field).name
        if not isinstance(self.alpha, Number):
            self.alpha = ds._get_field_info(self.alpha).name
            if self.luminosity_field[0] != self.alpha[0]:
                mylog.warning(
                    "The 'luminosity_field' %s and the 'alpha' field %s do not have the same field type!",
                    self.luminosity_field,
                    self.alpha,
                )
        self.ftype = self.luminosity_field[0]
        if mode == "spectrum":
            self.setup_pbar(data_source, self.luminosity_field)

    def __repr__(self):
        rets = [
            "PowerLawSourceModel(\n",
            f"    e0={self.e0}\n",
            f"    emin={self.emin}\n",
            f"    emax={self.emax}\n",
            f"    luminosity_field={self.luminosity_field}\n",
            f"    alpha={self.alpha}\n",
            ")",
        ]
        return "".join(rets)

    def cleanup_model(self, mode):
        if mode == "spectrum":
            self.pbar.close()

    def make_spectrum(
        self,
        data_source,
        emin,
        emax,
        nbins,
        redshift=0.0,
        dist=None,
        cosmology=None,
        normal=None,
    ):
        """
        Using all the data in a yt data container, make a count rate spectrum in the source frame,
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
        normal : integer, string, or array-like, optional
            This is a line-of-sight direction along which the spectrum will be
            Doppler shifted using the velocity field in the object. This is
            only an option if the spectrum is calculated in the observer frame.
            Options are one of "x", "y", "z", 0, 1, 2, or an 3-element array-like
            object of floats to specify an off-axis normal vector.

        Returns
        -------
        :class:`~soxs.spectra.CountRateSpectrum` or :class:`~soxs.spectra.Spectrum`,
        depending on how the method is invoked.
        """
        normal = sanitize_normal(normal)
        if normal is not None:
            if redshift == 0.0 and dist is None:
                raise RuntimeError(
                    "Cannot use a normal vector for the line-of-sight "
                    "when a redshift or the distance is not specified! not specified or the distance "
                )
            data_source.set_field_parameter("axis", normal)
        shifting = normal is not None
        ebins = np.linspace(emin, emax, nbins + 1)
        spec = np.zeros(nbins)
        spectral_norm = 1.0
        self.setup_model("spectrum", data_source, redshift)
        for chunk in data_source.chunks([], "io"):
            chunk_data = self.process_data(
                "spectrum", chunk, spectral_norm, shifting=shifting, ebins=ebins
            )
            if chunk_data is not None:
                spec += chunk_data
        self.cleanup_model("spectrum")
        return self._make_spectrum(
            data_source.ds, ebins, spec, redshift, dist, cosmology
        )

    def process_data(
        self,
        mode,
        chunk,
        spectral_norm,
        ebins=None,
        emin=None,
        emax=None,
        shifting=False,
    ):

        num_cells = check_num_cells(self.ftype, chunk)

        if num_cells == 0:
            if mode in ["photons", "spectrum"]:
                return
            else:
                return np.array([])

        if mode in ["spectrum", "intensity", "photon_intensity"] and shifting:
            shift = self.compute_shift(chunk)
        else:
            shift = np.ones_like(chunk[self.luminosity_field].d)

        if isinstance(self.alpha, float):
            alpha = self.alpha * np.ones_like(chunk[self.luminosity_field].d)
        else:
            alpha = chunk[self.alpha].d

        if emin is not None and emax is not None:
            ei = emin
            ef = emax
        else:
            ei = self.emin.v
            ef = self.emax.v

        etoalpha = self.e0.v**alpha
        K_fac = self.emax.v ** (2.0 - alpha) - self.emin.v ** (2.0 - alpha)
        K_fac[alpha == 2] = np.log(self.emax.v / self.emin.v)
        K_fac *= etoalpha
        if np.any(alpha != 2):
            K_fac[alpha != 2] /= 2.0 - alpha[alpha != 2]

        try:
            K = chunk[self.luminosity_field].to_value("keV/s") / K_fac
        except UnitConversionError:
            raise ValueError('The "luminosity_field" must be in units of power!')

        if mode in ["photons", "photon_rate", "photon_intensity"]:

            Nph = (ef / shift) ** (1.0 - alpha) - (ei / shift) ** (1.0 - alpha)
            Nph[alpha == 1] = np.log(ef / ei)
            Nph *= K * etoalpha
            if np.any(alpha != 1):
                Nph[alpha != 1] /= 1.0 - alpha[alpha != 1]

            if mode == "photons":
                Nph *= spectral_norm * self.scale_factor
                if self.observer == "internal":
                    r2 = self.compute_radius(chunk)
                    Nph /= r2

                number_of_photons = self.prng.poisson(lam=Nph)

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
                            e = ei ** (1.0 - alpha[i]) + u * Nph[i]
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

            else:
                return Nph * shift * shift

        elif mode in ["luminosity", "intensity"]:
            L = (ef / shift) ** (2.0 - alpha) - (ei / shift) ** (2.0 - alpha)
            L[alpha == 2] = np.log(ef / ei)
            L *= K * etoalpha
            if np.any(alpha != 2):
                L[alpha != 2] /= 2.0 - alpha[alpha != 2]

            return L * shift * shift * shift

        elif mode == "spectrum":
            inv_sf = 1.0 / self.scale_factor
            emid = 0.5 * (ebins[1:] + ebins[:-1]) * inv_sf / self.e0.v

            spec = power_law_spectrum(num_cells, emid, alpha, K, shift, self.pbar)

            return spec
