from numbers import Number

import numpy as np
from more_itertools import chunked
from soxs.constants import atomic_weights, elem_names, metal_elem
from soxs.utils import parse_prng
from unyt.array import unyt_quantity
from yt.data_objects.static_output import Dataset
from yt.utilities.exceptions import YTFieldNotFound

from pyxsim.lib.spectra import make_band, shift_spectrum
from pyxsim.source_models.data_handlers import find_data_handler
from pyxsim.source_models.sources import SourceModel
from pyxsim.utils import (
    _parse_abund_table,
    compute_H_abund,
    isunitful,
    mylog,
    parse_value,
)

keV_per_K = unyt_quantity(1.0, "K").to_value("keV", "thermal")


class ThermalSourceModel(SourceModel):
    _density_dependence = False
    _nei = False

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
        self.he_d_fraction = None  # Will be set by the subclass
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
        self._efluxf = None
        self._pfluxf = None

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
        self._efluxf = None
        self._pfluxf = None
        self.data_handler = find_data_handler(data_source)
        if isinstance(data_source, Dataset):
            ds = data_source
        else:
            ds = data_source.ds
        try:
            err_msg = f"The {self.emission_measure_field} field is not "
            "found, probably because the individual fields "
            "for hydrogen nuclei density and electron number "
            "density are not present. If you do not have species "
            "fields in your dataset, you may need to set "
            "default_species_fields='ionized' in the call "
            "to yt.load(), set them up using Trident, or "
            "set the field manually."
            self.emission_measure_field = self.data_handler.get_field_info(self.emission_measure_field)
            ftype = self.emission_measure_field.name[0]
        except YTFieldNotFound as e:
            raise RuntimeError(err_msg) from e
        self.temperature_field = self.data_handler.get_field_info(self.temperature_field)
        fields = [self.emission_measure_field.name, self.temperature_field.name]
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
        if self.he_d_fraction is not None:
            fields.append(self.he_d_fraction)
        ftypes = np.array([f[0] for f in fields])
        if not np.all(ftypes == ftype):
            mylog.warning("Not all fields have the same field type! Fields used: %s", fields)
        self.density_field = (ftype, "density")
        self.entropy_field = (ftype, "entropy")
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
        if self.he_d_fraction is not None:
            mylog.info("Using he_d_fraction '%s'.", self.he_d_fraction)
        self.spectral_model.prepare_spectrum(redshift)
        if mode in ["photons", "spectrum"]:
            self.setup_pbar(data_source, self.temperature_field)

    def fluxf(self, mode):
        if mode in ["luminosity", "intensity"]:
            if self._efluxf is None:
                self._efluxf = self.spectral_model.make_fluxf(self.emin, self.emax, energy=True)
            return self._efluxf
        elif mode in ["photon_rate", "photon_intensity"]:
            if self._pfluxf is None:
                self._pfluxf = self.spectral_model.make_fluxf(self.emin, self.emax, energy=False)
            return self._pfluxf
        else:
            raise NotImplementedError

    def _process_chunk(self, chunk, mode, shifting):
        out_chunk = {
            "orig_shape": chunk[self.density_field].shape,
            "density": np.ravel(self.data_handler.process_array(chunk[self.density_field], "g/cm**3")),
            "kT": np.ravel(keV_per_K * self.data_handler.process_array(chunk[self.temperature_field], "K")),
            "entropy": np.ravel(self.data_handler.process_array(chunk[self.entropy_field], "keV*cm**2")),
            "emission_measure": np.ravel(
                self.data_handler.process_array(chunk[self.emission_measure_field], "cm**-3")
            ),
        }
        num_cells = out_chunk["density"].size
        if mode in ["spectrum", "intensity", "photon_intensity"] and shifting:
            out_chunk["velocity_magnitude"] = self.data_handler.process_array(
                chunk[self.ftype, "velocity_magnitude"], "c"
            )
            out_chunk["velocity_los"] = self.data_handler.process_array(
                chunk[self.ftype, "velocity_los"], "c"
            )
        if self.nh_field is not None:
            out_chunk["H_nuclei_density"] = np.ravel(
                self.data_handler.process_array(chunk[self.nh_field], "1/cm**3")
            )
        if isinstance(self.h_fraction, Number):
            X_H = self.h_fraction
        else:
            X_H = np.ravel(self.data_handler.process_array(chunk[self.h_fraction]))
        if self._nei:
            out_chunk["metallicity"] = np.zeros(num_cells)
            elem_keys = self.var_ion_keys
        else:
            elem_keys = self.var_elem_keys
            if isinstance(self.Zmet, Number):
                out_chunk["metallicity"] = self.Zmet * np.ones(num_cells)
            else:
                out_chunk["metallicity"] = np.ravel(self.data_handler.process_array(chunk[self.Zmet]))
                fac = self.Zconvert
                if str(chunk[self.Zmet].units) != "Zsun":
                    fac /= X_H
                out_chunk["metallicity"] *= fac

        if self.num_var_elem > 0:
            for key in elem_keys:
                value = self.var_elem[key]
                if isinstance(value, Number):
                    out_chunk[f"{key}_abundance"] = value * np.ones(num_cells)
                else:
                    eZ = np.ravel(self.data_handler.process_array(chunk[value]))
                    fac = self.mconvert[key]
                    if str(chunk[value].units) != "Zsun":
                        fac /= X_H
                    out_chunk[f"{key}_abundance"] = eZ * fac
        return out_chunk

    def _process_data(
        self,
        mode,
        chunk,
        spectral_norm,
        ebins=None,
        emin=None,
        emax=None,
        shifting=False,
    ):
        if mode == "spectrum":
            spec = np.zeros(ebins.size - 1)
        else:
            spec = None

        shifted_intensity = mode.endswith("intensity") and shifting

        orig_shape = chunk["orig_shape"]
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
            cut &= chunk["density"] < self.max_density
        if self.min_entropy is not None:
            cut &= chunk["entropy"] > self.min_entropy
        kT = chunk["kT"]
        cut &= (kT >= self.kT_min) & (kT <= self.kT_max)
        metalZ = chunk["metallicity"]
        cell_nrm = chunk["emission_measure"] * spectral_norm

        nH = chunk.get("H_nuclei_density", None)

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
            return np.zeros(orig_shape)

        if mode in ["spectrum", "intensity", "photon_intensity"] and shifting:
            shift = self.compute_shift(chunk, cut=cut)
        else:
            shift = np.ones(num_cells)

        kT = kT[cut]
        cell_nrm = cell_nrm[cut]
        metalZ = metalZ[cut]
        elem_keys = self.var_ion_keys if self._nei else self.var_elem_keys
        if self.num_var_elem > 0:
            elemZ = np.zeros((self.num_var_elem, num_cells))
            for i, key in enumerate(elem_keys):
                elemZ[i] = chunk[f"{key}_abundance"][cut]
        if nH:
            nH = nH[cut]

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
            if self._density_dependence:
                nHi = nH[ibegin:iend]
            else:
                nHi = None

            if mode in ["photons", "spectrum"] or shifted_intensity:
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
                if self._density_dependence:
                    cflux, mflux, vflux = self.fluxf(mode)(kTi, nHi)
                else:
                    cflux, mflux, vflux = self.fluxf(mode)(kTi)
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
            return spec / np.diff(ebins)
        else:
            return np.resize(ret, orig_shape)
