from numbers import Number

import numpy as np
from more_itertools import chunked
from soxs.constants import atomic_weights, elem_names, metal_elem
from soxs.utils import parse_prng
from unyt.array import unyt_quantity
from yt.data_objects.static_output import Dataset
from yt.utilities.exceptions import YTFieldNotFound

from pyxsim.lib.spectra import make_band, shift_spectrum
from pyxsim.source_models.sources import SourceModel
from pyxsim.utils import (
    _parse_abund_table,
    compute_H_abund,
    isunitful,
    mylog,
    parse_value,
    sanitize_normal,
)


class ThermalSourceModel(SourceModel):
    _density_dependence = False
    _nei = False
    _cx = False

    def __init__(
        self,
        spectral_model,
        emin,
        emax,
        nbins,
        Zmet,
        binscale="linear",
        kT_min=0.025,
        kT_max=64.0,
        var_elem=None,
        max_density=None,
        min_entropy=None,
        method="invert_cdf",
        abund_table="angr",
        prng=None,
        temperature_field=None,
        emission_measure_field=None,
        h_fraction=None,
        nH_min=None,
        nH_max=None,
    ):
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
        self.trace_abund = None  # Will be set by the subclass
        if max_density is not None:
            if not isinstance(max_density, unyt_quantity):
                if isinstance(max_density, tuple):
                    max_density = unyt_quantity(max_density[0], max_density[1])
                else:
                    max_density = unyt_quantity(max_density, "g/cm**3")
        if min_entropy is not None:
            if not isinstance(min_entropy, unyt_quantity):
                if isinstance(min_entropy, tuple):
                    min_entropy = unyt_quantity(min_entropy[0], min_entropy[1])
                else:
                    min_entropy = unyt_quantity(min_entropy, "keV*cm**2")
        self.temperature_field = temperature_field
        self.emission_measure_field = emission_measure_field
        self.density_field = None  # Will be determined later
        self.nh_field = None  # Will be set by the subclass
        self.collnpar = None  # Will be set by the subclass
        self.h_r_number_density = None  # Will be set by the subclass
        self.h_d_number_density = None  # Will be set by the subclass
        self.he_d_number_density = None  # Will be set by the subclass
        self.max_density = max_density
        self.min_entropy = min_entropy
        self.tot_num_cells = 0  # Will be determined later
        self.ftype = "gas"
        self.binscale = binscale
        self.abund_table = abund_table
        self.method = method
        self.prng = parse_prng(prng)
        self.kT_min = kT_min
        self.kT_max = kT_max
        mylog.info("kT_min = %g keV", kT_min)
        mylog.info("kT_max = %g keV", kT_max)
        self.nH_min = nH_min
        self.nH_max = nH_max
        self.redshift = None
        self.pbar = None
        self.Zconvert = 1.0
        self.mconvert = {}
        self.abund_table = abund_table
        self.atable = _parse_abund_table(abund_table)
        if h_fraction is None:
            h_fraction = compute_H_abund(abund_table)
        self.h_fraction = h_fraction
        self.ebins = self.spectral_model.ebins
        self.de = self.spectral_model.de
        self.emid = self.spectral_model.emid
        self.bin_edges = np.log10(self.ebins) if self.binscale == "log" else self.ebins
        self.nbins = self.emid.size
        self.model_vers = self.spectral_model.model_vers
        self.efluxf = None
        self.pfluxf = None

    def _prep_repr(self):
        class_name = self.__class__.__name__
        strs = {
            "emin": self.emin,
            "emax": self.emax,
            "nbins": self.nbins,
            "Zmet": self.Zmet,
            "binscale": self.binscale,
            "temperature_field": self.temperature_field,
            "emission_measure_field": self.emission_measure_field,
            "kT_min": self.kT_min,
            "kT_max": self.kT_max,
            "method": self.method,
            "model_vers": self.spectral_model.model_vers,
            "max_density": self.max_density,
            "min_entropy": self.min_entropy,
            "abund_table": self.abund_table,
            "h_fraction": self.h_fraction,
            "var_elem": self.var_elem,
        }
        return class_name, strs

    def __repr__(self):
        class_name, strs = self._prep_repr()
        ret = f"{class_name}(\n"
        for key, value in strs.items():
            ret += f"    {key}={value}\n"
        ret += ")\n"
        return ret

    def setup_model(self, mode, data_source, redshift):
        if isinstance(data_source, Dataset):
            ds = data_source
        else:
            ds = data_source.ds
        try:
            if self._cx:
                err_msg = (
                    "One of the fields necessary to create the charge exchange "
                    "emision measure field is not present!"
                )
                self.h_r_number_density = ds._get_field_info(self.h_r_number_density).name
                self.h_d_number_density = ds._get_field_info(self.h_d_number_density).name
                self.he_d_number_density = ds._get_field_info(self.he_d_number_density).name
                ftype = self.h_r_number_density[0]

                def _emission_measure_cx(field, data):
                    dV = data[ftype, "mass"] / data[ftype, "density"]
                    n_h_r = data[self.h_r_number_density]
                    n_d = data[self.h_d_number_density] + data[self.he_d_number_density]
                    return n_h_r * n_d * dV

                ds.add_field(
                    (ftype, "emission_measure_cx"),
                    _emission_measure_cx,
                    units="cm**-3",
                    sampling_type="local",
                    force_override=True,
                )
                self.emission_measure_field = ds._get_field_info((ftype, "emission_measure_cx")).name
            else:
                err_msg = f"The {self.emission_measure_field} field is not "
                "found, probably because the individual fields "
                "for hydrogen nuclei density and electron number "
                "density are not present. If you do not have species "
                "fields in your dataset, you may need to set "
                "default_species_fields='ionized' in the call "
                "to yt.load(), set them up using Trident, or "
                "set the field manually."
                self.emission_measure_field = ds._get_field_info(self.emission_measure_field).name
                ftype = self.emission_measure_field[0]
        except YTFieldNotFound as e:
            raise RuntimeError(err_msg) from e
        self.temperature_field = ds._get_field_info(self.temperature_field).name
        fields = [self.emission_measure_field, self.temperature_field]
        self.ftype = ftype
        self.redshift = redshift
        if not self._nei and not isinstance(self.Zmet, Number):
            zfield = ds._get_field_info(self.Zmet)
            Z_units = str(zfield.units)
            self.Zmet = zfield.name
            fields.append(self.Zmet)
            if Z_units in ["dimensionless", "", "code_metallicity"]:
                Zsum = (self.atable * atomic_weights)[metal_elem].sum()
                self.Zconvert = atomic_weights[1] / Zsum
            elif Z_units == "Zsun":
                self.Zconvert = 1.0
            else:
                raise RuntimeError(f"I don't understand metallicity units of {Z_units}!")
        if self.num_var_elem > 0:
            for key in self.var_elem:
                value = self.var_elem[key]
                if not isinstance(value, Number):
                    if "^" in key:
                        elem = key.split("^")[0]
                    else:
                        elem = key
                    n_elem = elem_names.index(elem)
                    vfield = ds._get_field_info(value)
                    fields.append(vfield.name)
                    m_units = str(vfield.units)
                    self.var_elem[key] = vfield.name
                    if m_units in ["dimensionless", "", "code_metallicity"]:
                        m = self.atable[n_elem] * atomic_weights[n_elem]
                        self.mconvert[key] = atomic_weights[1] / m
                    elif m_units == "Zsun":
                        self.mconvert[key] = 1.0
                    else:
                        raise RuntimeError(f"I don't understand units of {m_units} for element {key}!")
        if self.nh_field is not None:
            self.nh_field = ds._get_field_info(self.nh_field).name
            fields.append(self.nh_field)
        if not isinstance(self.h_fraction, Number):
            self.h_fraction = ds._get_field_info(self.h_fraction).name
            fields.append(self.h_fraction)
        if self.h_r_number_density is not None:
            fields.append(self.h_r_number_density)
        if self.h_d_number_density is not None:
            fields.append(self.h_d_number_density)
        if self.he_d_number_density is not None:
            fields.append(self.h_d_number_density)
        ftypes = np.array([f[0] for f in fields])
        if not np.all(ftypes == ftype):
            mylog.warning("Not all fields have the same field type! Fields used: %s", fields)
        self.density_field = (ftype, "density")
        self.entropy_field = (ftype, "entropy")
        if not self._cx:
            mylog.info("Using emission measure field '%s'.", self.emission_measure_field)
        mylog.info("Using temperature field '%s'.", self.temperature_field)
        if self.nh_field is not None:
            mylog.info("Using nH field '%s'.", self.nh_field)
        if self.collnpar is not None:
            mylog.info("Using collnpar '%s'.", self.collnpar)
            if isunitful(self.collnpar):
                self.collnpar = float(parse_value(self.collnpar, "km/s").v)
        if self.h_r_number_density is not None:
            mylog.info("Using h_r_number_density '%s'.", self.h_r_number_density)
        if self.h_d_number_density is not None:
            mylog.info("Using h_d_number_density '%s'.", self.h_d_number_density)
        if self.he_d_number_density is not None:
            mylog.info("Using he_d_number_density '%s'.", self.he_d_number_density)
        self.spectral_model.prepare_spectrum(redshift)
        if mode in ["photons", "spectrum"]:
            self.setup_pbar(data_source, self.temperature_field)

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
        emin : float, (value, unit) tuple, unyt_quantity, or Quantity
            The minimum energy in the band. If a float, it is assumed to be
            in keV.
        emax : float, (value, unit) tuple, unyt_quantity, or Quantity
            The minimum energy in the band. If a float, it is assumed to be
            in keV.
        nbins : integer
            The number of bins in the spectrum.
        redshift : float, optional
            If greater than 0, we assume that the spectrum should be created in
            the observer frame at a distance given by the cosmology. Default: 0.0
        dist : float, (value, unit) tuple, unyt_quantity, or Quantity, optional
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
        self.setup_model("spectrum", data_source, redshift)
        spectral_norm = 1.0
        spec = np.zeros(nbins)
        ebins = np.linspace(emin, emax, nbins + 1)
        for chunk in data_source.chunks([], "io"):
            chunk_data = self.process_data("spectrum", chunk, spectral_norm, shifting=shifting, ebins=ebins)
            if chunk_data is not None:
                spec += chunk_data
        spec /= np.diff(ebins)
        self.cleanup_model("spectrum")
        return self._make_spectrum(data_source.ds, ebins, spec, redshift, dist, cosmology)

    def make_fluxf(self, emin, emax, energy=False):
        return self.spectral_model.make_fluxf(emin, emax, energy=energy)

    def process_data(
        self,
        mode,
        chunk,
        spectral_norm,
        ebins=None,
        emin=None,
        emax=None,
        fluxf=None,
        shifting=False,
    ):
        if mode == "spectrum":
            spec = np.zeros(ebins.size - 1)
        else:
            spec = None

        shifted_intensity = mode.endswith("intensity") and shifting

        orig_shape = chunk[self.temperature_field].shape
        if len(orig_shape) == 0:
            orig_ncells = 0
        else:
            orig_ncells = np.prod(orig_shape)
        if orig_ncells == 0:
            if mode in ["photons", "spectrum"]:
                return
            else:
                return np.array([])

        ret = np.zeros(orig_ncells)

        cut = True

        if self.max_density is not None:
            cut &= np.ravel(chunk[self.density_field]) < self.max_density
        if self.min_entropy is not None:
            cut &= np.ravel(chunk[self.entropy_field]) > self.min_entropy
        kT = np.ravel(chunk[self.temperature_field].to_value("keV", "thermal"))
        cut &= (kT >= self.kT_min) & (kT <= self.kT_max)

        cell_nrm = np.ravel(chunk[self.emission_measure_field].d * spectral_norm)

        if self.nh_field is not None:
            nH = np.ravel(chunk[self.nh_field].d)
        else:
            nH = None

        if self._cx:
            if isinstance(self.collnpar, Number):
                coll = self.collnpar
            else:
                coll = np.ravel(chunk[self.collnpar].d)
            n_h_d = np.ravel(chunk[self.h_d_number_density].d)
            n_he_d = np.ravel(chunk[self.he_d_number_density].d)
        else:
            coll = None
            n_h_d = None
            n_he_d = None

        if isinstance(self.h_fraction, Number):
            X_H = self.h_fraction
        else:
            X_H = np.ravel(chunk[self.h_fraction].d)

        num_cells = cut.sum()

        if mode in ["photons", "spectrum"]:
            if num_cells == 0:
                self.pbar.update(orig_ncells)
                return
            else:
                self.pbar.update(orig_ncells - num_cells)
        elif num_cells == 0:
            # Here, we have no active cells, and so we
            # return an array of zeros with the original shape.
            # But yt needs to know that we may depend on various
            # fields, so we check for them here. Very hacky!
            if not isinstance(self.Zmet, Number):
                _ = chunk[self.Zmet]
            if self._cx:
                if not isinstance(self.collnpar, Number):
                    _ = chunk[self.collnpar]
            if self.num_var_elem > 0:
                elem_keys = self.var_ion_keys if self._nei else self.var_elem_keys
                for key in elem_keys:
                    value = self.var_elem[key]
                    if not isinstance(value, Number):
                        _ = chunk[value]
            # We also need to do this for the velocity fields if we use them
            if mode in ["spectrum", "intensity", "photon_intensity"] and shifting:
                _ = chunk[self.ftype, "velocity_magnitude"]
            return np.zeros(orig_shape)

        if mode in ["spectrum", "intensity", "photon_intensity"] and shifting:
            shift = self.compute_shift(chunk, cut=cut)
        else:
            shift = np.ones(num_cells)

        kT = kT[cut]
        cell_nrm = cell_nrm[cut]
        if nH is not None:
            nH = nH[cut]
        if self._cx:
            if isinstance(self.collnpar, Number):
                coll = coll * np.ones(num_cells)
            else:
                coll = coll[cut]
            n_h_d = n_h_d[cut]
            n_he_d = n_he_d[cut]

        if not isinstance(X_H, Number):
            X_H = X_H[cut]

        if self._nei:
            metalZ = np.zeros(num_cells)
            elem_keys = self.var_ion_keys
        else:
            elem_keys = self.var_elem_keys
            if isinstance(self.Zmet, Number):
                metalZ = self.Zmet * np.ones(num_cells)
            else:
                mZ = chunk[self.Zmet]
                fac = self.Zconvert
                if str(mZ.units) != "Zsun":
                    fac /= X_H
                metalZ = np.ravel(mZ.d * fac)[cut]

        elemZ = None
        if self.num_var_elem > 0:
            elemZ = np.zeros((self.num_var_elem, num_cells))
            for j, key in enumerate(elem_keys):
                value = self.var_elem[key]
                if isinstance(value, Number):
                    elemZ[j, :] = value
                else:
                    eZ = chunk[value]
                    fac = self.mconvert[key]
                    if str(eZ.units) != "Zsun":
                        fac /= X_H
                    elemZ[j, :] = np.ravel(eZ.d * fac)[cut]

        if self.observer == "internal" and mode == "photons":
            r2 = self.compute_radius(chunk, cut=cut)
            cell_nrm /= r2

        num_photons_max = 10000000
        number_of_photons = np.zeros(num_cells, dtype="int64")
        energies = np.zeros(num_photons_max)

        start_e = 0
        end_e = 0

        idxs = np.where(cut)[0]

        for ck in chunked(range(num_cells), 100):
            ibegin = ck[0]
            iend = ck[-1] + 1
            nck = iend - ibegin

            cnm = cell_nrm[ibegin:iend]

            kTi = kT[ibegin:iend]

            shifti = shift[ibegin:iend]
            if self._cx:
                colli = coll[ibegin:iend]
                with np.errstate(divide="ignore", invalid="ignore"):
                    h_f = n_h_d[ibegin:iend] / (n_h_d[ibegin:iend] + n_he_d[ibegin:iend])
                    he_f = 1.0 - h_f
                h_f = np.nan_to_num(h_f)
                he_f = np.nan_to_num(he_f)
            else:
                colli = None
                h_f = None
                he_f = None
            if self._density_dependence:
                nHi = nH[ibegin:iend]
            else:
                nHi = None

            if mode in ["photons", "spectrum"] or shifted_intensity:
                if self._cx:
                    tot_spec = np.zeros((self.nck, self.nbins))
                    if self._nei:
                        h_mspec, he_mspec, h_vspec, he_vspec = self.spectral_model.get_spectrum(colli)
                    else:
                        h_mspec, he_mspec, h_vspec, he_vspec = self.spectral_model.get_spectrum(kTi, colli)
                        tot_spec += metalZ[ibegin:iend, np.newaxis] * h_f * h_mspec + he_f * he_mspec
                    if self.num_var_elem > 0:
                        tot_spec += np.sum(
                            elemZ[:, ibegin:iend, np.newaxis] * h_f * h_vspec + he_f * he_vspec,
                            axis=0,
                        )
                else:
                    if self._density_dependence:
                        cspec, mspec, vspec = self.spectral_model.get_spectrum(kTi, nHi)
                    else:
                        cspec, mspec, vspec = self.spectral_model.get_spectrum(kTi)
                    tot_spec = cspec
                    tot_spec += metalZ[ibegin:iend, np.newaxis] * mspec
                    if self.num_var_elem > 0:
                        tot_spec += np.sum(elemZ[:, ibegin:iend, np.newaxis] * vspec, axis=0)
                np.clip(tot_spec, 0.0, None, out=tot_spec)

                if mode == "photons":
                    spec_sum = tot_spec.sum(axis=-1)
                    cell_norm = spec_sum * cnm

                    cell_n = np.atleast_1d(self.prng.poisson(lam=cell_norm))

                    number_of_photons[ibegin:iend] = cell_n
                    end_e += int(cell_n.sum())

                    norm_factor = 1.0 / spec_sum
                    p = norm_factor[:, np.newaxis] * tot_spec
                    cp = np.insert(np.cumsum(p, axis=-1), 0, 0.0, axis=1)
                    ei = start_e
                    for icell in range(nck):
                        cn = cell_n[icell]
                        if cn == 0:
                            continue
                        if self.method == "invert_cdf":
                            randvec = self.prng.uniform(size=cn)
                            randvec.sort()
                            cell_e = np.interp(randvec, cp[icell, :], self.bin_edges)
                        elif self.method == "accept_reject":
                            eidxs = self.prng.choice(self.nbins, size=cn, p=p[icell, :])
                            cell_e = self.emid[eidxs]
                        while ei + cn > num_photons_max:
                            num_photons_max *= 2
                        if num_photons_max > energies.size:
                            energies.resize(num_photons_max, refcheck=False)
                        energies[ei : ei + cn] = cell_e
                        ei += cn
                    start_e = end_e

                elif mode == "spectrum":
                    spec += shift_spectrum(self.ebins, ebins, tot_spec, shifti, cnm)

                elif mode.endswith("intensity"):
                    use_energy = int(mode == "intensity")
                    I = make_band(use_energy, emin, emax, self.ebins, self.emid, tot_spec, shift)
                    ret[idxs[ibegin:iend]] = I * cnm

                if mode in ["photons", "spectrum"]:
                    self.pbar.update(nck)

            else:
                if self._cx:
                    if self._nei:
                        h_flux, he_flux = fluxf(colli)
                    else:
                        h_flux, he_flux = fluxf(kT, colli)
                    hhe_flux = h_f * h_flux + he_f * he_flux
                    tot_flux = metalZ[ibegin:iend] * hhe_flux
                    if self.num_var_elem > 0:
                        tot_flux += np.sum(elemZ[:, ibegin:iend] * hhe_flux, axis=0)
                else:
                    if self._density_dependence:
                        cflux, mflux, vflux = fluxf(kTi, nHi)
                    else:
                        cflux, mflux, vflux = fluxf(kTi)
                    tot_flux = cflux
                    tot_flux += metalZ[ibegin:iend] * mflux
                    if self.num_var_elem > 0:
                        tot_flux += np.sum(elemZ[:, ibegin:iend] * vflux, axis=0)

                ret[idxs[ibegin:iend]] = tot_flux * cnm

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

    def cleanup_model(self, mode):
        if mode in ["spectrum", "photons"]:
            self.pbar.close()
