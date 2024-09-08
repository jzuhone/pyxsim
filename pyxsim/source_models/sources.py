import numpy as np
from soxs.spectra import CountRateSpectrum, Spectrum
from soxs.utils import parse_prng
from tqdm.auto import tqdm
from unyt.array import unyt_quantity
from yt.utilities.cosmology import Cosmology
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    communication_system,
    parallel_capable,
    parallel_objects,
)

from pyxsim.utils import ParallelProgressBar, parse_value

cm2_per_kpc2 = unyt_quantity(1.0, "kpc**2").to_value("cm**2")

comm = communication_system.communicators[-1]


class SourceModel:
    def __init__(self, prng=None):
        self.spectral_norm = None
        self.redshift = None
        self.prng = parse_prng(prng)
        self.observer = "external"

    def process_data(self, mode, chunk, spectral_norm, fluxf=None):
        # This needs to be implemented for every
        # source model specifically
        pass

    def setup_pbar(self, data_source, field):
        citer = data_source.chunks([], "io")
        num_cells = 0
        for chunk in parallel_objects(citer):
            num_cells += chunk[field].size
        self.tot_num_cells = comm.mpi_allreduce(num_cells)
        if parallel_capable:
            self.pbar = ParallelProgressBar("Processing cells/particles ")
        else:
            self.pbar = tqdm(
                leave=True, total=self.tot_num_cells, desc="Processing cells/particles "
            )

    def setup_model(self, mode, data_source, redshift):
        # This needs to be implemented for every
        # source model specifically
        pass

    def set_pv(self, p_fields, v_fields, le, re, dw, c, periodicity, observer):
        self.p_fields = p_fields
        self.v_fields = v_fields
        self.le = le
        self.re = re
        self.dw = dw
        self.c = c
        self.periodicity = periodicity
        self.observer = observer

    def compute_radius(self, chunk, cut=None):
        if cut is None:
            cut = ...
        pos = np.array(
            [np.ravel(chunk[self.p_fields[i]].to_value("kpc"))[cut] for i in range(3)]
        )
        for i in range(3):
            if self.periodicity[i]:
                tfl = pos[i] < self.le[i]
                tfr = pos[i] > self.re[i]
                pos[:, tfl] += self.dw[i]
                pos[:, tfr] -= self.dw[i]
        return np.sum((pos - self.c[:, np.newaxis]) ** 2, axis=0) * cm2_per_kpc2

    def cleanup_model(self, mode):
        # This needs to be implemented for every
        # source model specifically
        pass

    def make_fluxf(self, emin, emax, energy=False):
        # This needs to be implemented for every
        # source model specifically
        pass

    def _make_dist_fac(self, ds, redshift, dist, cosmology, per_sa=False):
        if dist is None:
            if cosmology is None:
                if hasattr(ds, "cosmology"):
                    cosmology = ds.cosmology
                else:
                    cosmology = Cosmology()
            dist = cosmology.luminosity_distance(0.0, redshift)
            angular_scale = 1.0 / cosmology.angular_scale(0.0, redshift)
        else:
            redshift = 0.0  # Only for local sources!
            try:
                # normal behaviour, if dist is a YTQuantity
                dist = ds.quan(dist.value, dist.units)
            except AttributeError as e:
                try:
                    dist = ds.quan(*dist)
                except (RuntimeError, TypeError):
                    raise TypeError(
                        "dist should be a YTQuantity or a (value, unit) tuple!"
                    ) from e
            angular_scale = dist / ds.quan(1.0, "radian")
        dist_fac = 1.0 / (4.0 * np.pi * dist * dist)
        if per_sa:
            dist_fac = ds.quan((dist_fac / angular_scale**2).v, "rad**-2")
        return dist_fac, redshift

    def _make_spectrum(self, ds, ebins, spec, redshift, dist, cosmology):
        if redshift > 0.0 or dist is not None:
            dist_fac, redshift = self._make_dist_fac(ds, redshift, dist, cosmology)
            spec *= (1.0 + redshift) ** 2 * dist_fac.in_cgs()
            spec_class = Spectrum
        else:
            spec_class = CountRateSpectrum
        return spec_class(ebins, spec)

    def make_source_fields(self, ds, emin, emax, force_override=False, band_name=None):
        """
        Make the following fields in the rest frame of the
        source within a specific energy band for a dataset in yt:

        f"xray_emissivity_{emin}_{emax}_keV" (in erg/cm**3/s)
        f"xray_luminosity_{emin}_{emax}_keV" (in erg/s)
        f"xray_photon_emissivity_{emin}_{emax}_keV" (in photons/cm**3/s)
        f"xray_count_rate_{emin}_{emax}_keV" (in photons/s)

        where "emin" and "emax" are the bounds of the energy band
        as described below. In this case, emin and emax are the bounds
        of the energy in the source frame.


        Parameters
        ----------
        ds : :class:`~yt.data_objects.static_output.Dataset`
            The loaded yt dataset to make the fields for.
        emin : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The minimum energy in the band in the source frame.
            If a float, it is assumed to be in keV.
        emax : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The minimum energy in the band in the source frame.
            If a float, it is assumed to be in keV.
        force_override : boolean, optional
            If True, override a pre-existing field with the same name.
            Default: False
        band_name : string, optional
            The name to give to the energy band in the field. If None,
            it is set to "{emin}_{emax}_keV". Default: None

        Returns
        -------
        The list of fields which are generated.
        """
        spectral_norm = 1.0
        redshift = 0.0

        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")

        self.setup_model("fields", ds, redshift)

        ftype = self.ftype

        if band_name is None:
            band_name = f"{emin.value}_{emax.value}_keV"

        emiss_name = (ftype, f"xray_emissivity_{band_name}")
        emiss_dname = rf"\epsilon_{{X}} ({emin.value:.2f}-{emax.value:.2f} keV)"

        lum_name = (ftype, emiss_name[1].replace("emissivity", "luminosity"))
        lum_dname = emiss_dname.replace(r"\epsilon", r"\rm{{L}}")

        phot_emiss_name = (
            ftype,
            emiss_name[1].replace("emissivity", "photon_emissivity"),
        )

        count_rate_name = (
            ftype,
            phot_emiss_name[1].replace("emissivity", "count_rate"),
        )
        count_rate_dname = lum_dname.replace("\rm{{L}}", "\rm{{R}}")

        efluxf = self.make_fluxf(emin, emax, energy=True)

        def _luminosity_field(field, data):
            return data.ds.arr(
                self.process_data("energy_field", data, spectral_norm, fluxf=efluxf),
                "keV/s",
            )

        ds.add_field(
            lum_name,
            function=_luminosity_field,
            display_name=lum_dname,
            sampling_type="local",
            units="erg/s",
            force_override=force_override,
        )

        def _emissivity_field(field, data):
            ret = data[lum_name]
            return ret * data[ftype, "density"] / data[ftype, "mass"]

        ds.add_field(
            emiss_name,
            function=_emissivity_field,
            display_name=emiss_dname,
            sampling_type="local",
            units="erg/cm**3/s",
            force_override=force_override,
        )

        pfluxf = self.make_fluxf(emin, emax, energy=False)

        def _count_rate_field(field, data):
            return data.ds.arr(
                self.process_data("photon_field", data, spectral_norm, fluxf=pfluxf),
                "photons/s",
            )

        ds.add_field(
            count_rate_name,
            function=_count_rate_field,
            display_name=count_rate_dname,
            sampling_type="local",
            units="photons/s",
            force_override=force_override,
        )

        def _photon_emissivity_field(field, data):
            ret = data[count_rate_name]
            return ret * data[ftype, "density"] / data[ftype, "mass"]

        ds.add_field(
            phot_emiss_name,
            function=_photon_emissivity_field,
            display_name=emiss_dname,
            sampling_type="local",
            units="photons/cm**3/s",
            force_override=force_override,
        )

        return [emiss_name, lum_name, phot_emiss_name, count_rate_name]

    def make_line_source_fields(self, ds, e0, de, line_name, force_override=False):
        """
        Make a list of source fields for a very small bandpass,
        essentially covering a line. The fields are of the form:

        f"xray_emissivity_{line_name}" (in erg/cm**3/s)
        f"xray_luminosity_{line_name}" (in erg/s)
        f"xray_photon_emissivity_{line_name}" (in photons/cm**3/s)
        f"xray_count_rate_{emin}_{line_name}" (in photons/s)

        In this case, e0 and de are in the source frame.

        Parameters
        ----------
        ds : :class:`~yt.data_objects.static_output.Dataset`
            The loaded yt dataset to make the fields for.
        e0 : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The line centroid in the source frame. If a float, it is
            assumed to be in keV.
        de : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The width of the band in the source frame. If a float, it is
            assumed to be in keV.
        band_name : string
            The suffix of the field name for the line. Logical names are
            "O_VII", "O_VIII", "Ne_IX", etc.
        force_override : boolean, optional
            If True, override a pre-existing field with the same name.
            Default: False

        Returns
        -------
        The list of fields which are generated.
        """
        e0 = parse_value(e0, "keV")
        de = parse_value(de, "keV")

        emin = e0 - 0.5 * de
        emax = e0 + 0.5 * de

        return self.make_source_fields(
            ds, emin, emax, band_name=line_name, force_override=force_override
        )

    def make_intensity_fields(
        self,
        ds,
        emin,
        emax,
        redshift=0.0,
        dist=None,
        cosmology=None,
        force_override=True,
        band_name=None,
    ):
        """
        Make the following fields in the observer frame within a
        specific energy band for a dataset in yt:

        f"xray_intensity_{emin}_{emax}_keV" (in erg/cm**3/s/arcsec**2)
        f"xray_photon_intensity_{emin}_{emax}_keV" (in photons/cm**3/s/arcsec**2)

        where "emin" and "emax" are the bounds of the energy band as
        described below. These should mainly be used for projections.
        In this case, emin and emax are the bounds of the energy in the
        observer frame.

        Parameters
        ----------
        ds : :class:`~yt.data_objects.static_output.Dataset`
            The loaded yt dataset to make the fields for.
        emin : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The minimum energy in the band in the observer frame. If a float,
            it is assumed to be in keV.
        emax : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The maximum energy in the band in the observer frame. If a float,
            it is assumed to be in keV.
        redshift : float, optional
            The redshift of the source. Default: 0.0
        dist : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The angular diameter distance, used for nearby sources. This may be
            optionally supplied instead of it being determined from the
            *redshift* and given *cosmology*. If units are not specified, it is
            assumed to be in kpc. To use this, the redshift must be set to zero.
        cosmology : :class:`~yt.utilities.cosmology.Cosmology`, optional
            Cosmological information. If not supplied, we try to get the
            cosmology from the dataset. Otherwise, LCDM with the default yt
            parameters is assumed.
        force_override : boolean, optional
            If True, override a pre-existing field with the same name.
            Default: False
        band_name : string, optional
            The name to give to the energy band in the field. If None,
            it is set to "{emin}_{emax}_keV". Default: None

        Returns
        -------
        The list of fields which are generated.
        """
        if redshift == 0.0 and dist is None:
            raise ValueError(
                "Either 'redshift' must be > 0.0 or 'dist' must " "not be None!"
            )

        spectral_norm = 1.0

        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")

        self.setup_model("fields", ds, 0.0)

        ftype = self.ftype

        dist_fac, redshift = self._make_dist_fac(
            ds, redshift, dist, cosmology, per_sa=True
        )

        emin_src = emin * (1.0 + redshift)
        emax_src = emax * (1.0 + redshift)

        if band_name is None:
            band_name = f"{emin.value}_{emax.value}_keV"

        ei_name = (ftype, f"xray_intensity_{band_name}")
        ei_dname = rf"I_{{X}} ({emin.value:.2f}-{emax.value:.2f} keV)"

        eif = self.make_fluxf(emin_src, emax_src, energy=True)

        def _intensity_field(field, data):
            ret = data.ds.arr(
                self.process_data("energy_field", data, spectral_norm, fluxf=eif),
                "keV/s",
            )
            idV = data[ftype, "density"] / data[ftype, "mass"]
            I = dist_fac * ret * idV
            return I.in_units("erg/cm**3/s/arcsec**2")

        ds.add_field(
            ei_name,
            function=_intensity_field,
            display_name=ei_dname,
            sampling_type="local",
            units="erg/cm**3/s/arcsec**2",
            force_override=force_override,
        )

        i_name = (ftype, ei_name[1].replace("intensity", "photon_intensity"))

        pif = self.make_fluxf(emin_src, emax_src, energy=False)

        def _photon_intensity_field(field, data):
            ret = data.ds.arr(
                self.process_data("photon_field", data, spectral_norm, fluxf=pif),
                "photons/s",
            )
            idV = data[ftype, "density"] / data[ftype, "mass"]
            I = (1.0 + redshift) * dist_fac * ret * idV
            return I.in_units("photons/cm**3/s/arcsec**2")

        ds.add_field(
            i_name,
            function=_photon_intensity_field,
            display_name=ei_dname,
            sampling_type="local",
            units="photons/cm**3/s/arcsec**2",
            force_override=force_override,
        )

        return [ei_name, i_name]

    def make_line_intensity_fields(
        self,
        ds,
        e0,
        de,
        line_name,
        redshift=0.0,
        dist=None,
        cosmology=None,
        force_override=False,
    ):
        """
        Make a list of intensity fields for a very small bandpass,
        essentially covering a line. The fields are of the form:

        f"xray_intensity_{line_name}" (in erg/cm**3/s/arcsec**2)
        f"xray_photon_intensity_{line_name}" (in photons/cm**3/s/arcsec**2)

        These should mainly be used for projections. In this case, e0 and de
        are in the observer frame.

        Parameters
        ----------
        ds : :class:`~yt.data_objects.static_output.Dataset`
            The loaded yt dataset to make the fields for.
        e0 : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The line centroid in the observer frame. If a float, it is
            assumed to be in keV.
        de : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The width of the band in the observer frame. If a float,
            it is assumed to be in keV.
        band_name : string
            The suffix of the field name for the line. Logical names are
            "O_VII", "O_VIII", "Ne_IX", etc.
        redshift : float, optional
            The redshift of the source. Default: 0.0
        dist : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The angular diameter distance, used for nearby sources. This may be
            optionally supplied instead of it being determined from the
            *redshift* and given *cosmology*. If units are not specified, it is
            assumed to be in kpc. To use this, the redshift must be set to zero.
        cosmology : :class:`~yt.utilities.cosmology.Cosmology`, optional
            Cosmological information. If not supplied, we try to get the
            cosmology from the dataset. Otherwise, LCDM with the default yt
            parameters is assumed.
        force_override : boolean, optional
            If True, override a pre-existing field with the same name.
            Default: False

        Returns
        -------
        The list of fields which are generated.
        """
        e0 = parse_value(e0, "keV")
        de = parse_value(de, "keV")

        emin = e0 - 0.5 * de
        emax = e0 + 0.5 * de

        return self.make_intensity_fields(
            ds,
            emin,
            emax,
            redshift=redshift,
            dist=dist,
            cosmology=cosmology,
            band_name=line_name,
            force_override=force_override,
        )
