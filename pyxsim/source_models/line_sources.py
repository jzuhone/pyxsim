from numbers import Number

import numpy as np
from scipy.stats import norm
from soxs.utils import parse_prng
from unyt.array import unyt_quantity
from yt.data_objects.static_output import Dataset
from yt.utilities.physical_constants import clight

from pyxsim.lib.spectra import line_spectrum
from pyxsim.source_models.sources import SourceModel
from pyxsim.utils import check_num_cells, isunitful, mylog, parse_value

gx = np.linspace(-7, 7, 10000)
gcdf = norm.cdf(gx)
gpdf = norm.pdf(gx)


class LineSourceModel(SourceModel):
    r"""
    Initialize a source model from a single line.

    Parameters
    ----------
    e0 : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The location of the emission line in energy in the rest frame of the
        source. If units are not given, they are assumed to be in keV.
    emission_field : string or (ftype, fname) tuple
        The field corresponding to the photon count rate per cell or particle,
        in the rest frame of the source, which serves as the normalization for
        the line. Must be in counts/s.
    sigma : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The standard intrinsic deviation of the emission line (not from Doppler
        broadening, which is handled in the projection step). Units of
        velocity or energy are accepted. If units are not given, they
        are assumed to be in keV. If set to a field name, the line broadening
        is assumed to be based on this field (in units of velocity or energy).
    prng : integer or :class:`~numpy.random.RandomState` object
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is to use the :mod:`numpy.random` module.

    Examples
    --------
    >>> location = (3.5, "keV")
    >>> sigma = (1000., "km/s")
    >>> line_model = LineSourceModel(location, "dark_matter_density_squared", sigma=sigma)
    """

    def __init__(self, e0, emission_field, sigma, prng=None):
        from unyt.exceptions import UnitConversionError

        self.e0 = parse_value(e0, "keV")
        if isinstance(sigma, Number):
            self.sigma = parse_value(sigma, "keV")
        elif isunitful(sigma):
            # The broadening is constant
            try:
                self.sigma = parse_value(sigma, "km/s")
                self.sigma *= self.e0 / clight
                self.sigma.convert_to_units("keV")
            except UnitConversionError:
                self.sigma = parse_value(sigma, "keV")
        else:
            # Should be a field name
            self.sigma = sigma
        self.emission_field = emission_field
        self.prng = parse_prng(prng)
        self.ftype = None

    def setup_model(self, mode, data_source, redshift):
        if isinstance(data_source, Dataset):
            ds = data_source
        else:
            ds = data_source.ds
        self.scale_factor = 1.0 / (1.0 + redshift)
        self.emission_field = ds._get_field_info(self.emission_field).name
        if not isinstance(self.sigma, (Number, unyt_quantity)):
            self.sigma = ds._get_field_info(self.sigma).name
            if self.emission_field[0] != self.sigma[0]:
                mylog.warning(
                    "The 'emission_field' %s and the 'sigma' field %s do not have the same field type!",
                    self.emission_field,
                    self.sigma,
                )
        self.ftype = self.emission_field[0]
        if mode == "spectrum":
            self.setup_pbar(data_source, self.emission_field)

    def __repr__(self):
        rets = [
            "LineSourceModel(\n",
            f"    e0={self.e0}\n",
            f"    emission_field={self.emission_field}\n",
            f"    sigma={self.sigma}\n",
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
        dist : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`, optional
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
            chunk_data = self.process_data(
                "spectrum", chunk, spectral_norm, ebins=ebins
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

        norm_field = chunk[self.emission_field]

        if mode in ["spectrum", "intensity", "photon_intensity"] and shifting:
            shift = self.compute_shift(chunk)
        else:
            shift = np.ones_like(norm_field.d)

        if isinstance(self.sigma, unyt_quantity):
            sigma = self.sigma.value * np.ones_like(norm_field.d)
        else:
            sigma = (chunk[self.sigma] * self.e0 / clight).to_value("keV")

        if mode == "photons":
            F = norm_field * spectral_norm * self.scale_factor
            if self.observer == "internal":
                r2 = self.compute_radius(chunk)
                F /= r2

            number_of_photons = self.prng.poisson(lam=F.in_cgs().d)

            energies = self.e0 * np.ones(number_of_photons.sum())

            start_e = 0
            for i in range(num_cells):
                if number_of_photons[i] > 0:
                    end_e = start_e + number_of_photons[i]
                    dE = (
                        self.prng.normal(
                            loc=0.0,
                            scale=sigma[i],
                            size=number_of_photons[i],
                        )
                        * self.e0.uq
                    )
                    energies[start_e:end_e] += dE
                    start_e = end_e

            energies = energies * self.scale_factor

            active_cells = number_of_photons > 0
            ncells = active_cells.sum()

            return ncells, number_of_photons[active_cells], active_cells, energies

        elif mode == "spectrum":

            ee = 0.5 * (ebins[1:] + ebins[:-1]) / self.scale_factor

            spec = line_spectrum(
                num_cells,
                float(self.e0),
                ee,
                sigma,
                gx,
                gpdf,
                norm_field.d,
                shift,
                self.pbar,
            )

            return spec

        else:

            xlo = emin - self.e0.value
            xhi = emax - self.e0.value
            xhis = xhi / sigma
            xlos = xlo / sigma
            fac = (norm.cdf(xhis) - norm.cdf(xlos)) * shift * shift
            if mode in ["luminosity", "intensity"]:
                fac = self.e0.value * fac
                fac -= (
                    sigma
                    * (np.exp(-0.5 * xhis**2) - np.exp(-0.5 * xlos**2))
                    / np.sqrt(2.0 * np.pi)
                )
                fac *= shift
            return fac * norm_field
