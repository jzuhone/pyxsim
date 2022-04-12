"""
Classes for specific source models
"""
import numpy as np
from yt.funcs import ensure_numpy_array
from tqdm.auto import tqdm
from pyxsim.utils import mylog
from yt.data_objects.static_output import Dataset
from yt.units.yt_array import YTQuantity, YTArray
from yt.utilities.physical_constants import clight
from yt.utilities.cosmology import Cosmology
from pyxsim.spectral_models import TableApecModel, XSpecAtableModel
from pyxsim.utils import parse_value, isunitful
from soxs.utils import parse_prng
from soxs.constants import elem_names, atomic_weights, metal_elem
from yt.utilities.exceptions import YTFieldNotFound
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects, communication_system, parallel_capable
from numbers import Number
from scipy.stats import norm
from more_itertools import chunked


gx = np.linspace(-6, 6, 2400)
gcdf = norm.cdf(gx)

comm = communication_system.communicators[-1]

solar_H_abund = 0.74

sqrt_two = np.sqrt(2.)


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


class SourceModel:
    def __init__(self, prng=None):
        self.spectral_norm = None
        self.redshift = None
        self.prng = parse_prng(prng)

    def process_data(self, mode, chunk):
        pass

    def setup_model(self, data_source, redshift, spectral_norm):
        self.spectral_norm = spectral_norm
        self.redshift = redshift

    def cleanup_model(self):
        pass

    def make_xray_fields(self, ds, redshift=0.0, dist=None, cosmology=None):
        r"""
        
        Parameters
        ----------
        ds : `~yt.data_objects.static_output.Dataset`
            The loaded yt dataset to make the fields for.
        redshift : float, optional
            The redshift of the source. Default: 0.0
        dist : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity, optional     
            The angular diameter distance, used for nearby sources. This may be
            optionally supplied instead of it being determined from the
            *redshift* and given *cosmology*. If units are not specified, it is
            assumed to be in kpc. To use this, the redshift must be set to zero.

        cosmology

        Returns
        -------
        The list of fields which are generated.
        """
        spectral_norm = 1.0

        self.setup_model(ds, redshift, spectral_norm)

        ftype = self.ftype

        def _emissivity_field(field, data):
            ret = data.ds.arr(self.process_data("energy_field", data), "keV/s")
            return ret*data[ftype, "density"]/data[ftype, "mass"]

        if hasattr(self, "emin"):
            emiss_name = (
                ftype, f"xray_emissivity_{self.emin.value}_{self.emax.value}_keV"
            )
            emiss_dname = rf"\epsilon_{{X}} ({self.emin.value}-{self.emax.value} keV)"
        else:
            emiss_name = (
                ftype, f"xray_emissivity_{self.e0.value}_keV"
            )
            emiss_dname = rf"\epsilon_{{X}} ({self.e0.value} keV)"

        ds.add_field(
            emiss_name,
            function=_emissivity_field,
            display_name=emiss_dname,
            sampling_type="local",
            units="erg/cm**3/s",
        )

        def _luminosity_field(field, data):
            return data.ds.arr(self.process_data("energy_field", data), "keV/s")

        lum_name = (
            ftype, emiss_name[1].replace("emissivity", "luminosity")
        )
        lum_dname = emiss_dname.replace("\epsilon", "\rm{{L}}")

        ds.add_field(
            lum_name,
            function=_luminosity_field,
            display_name=lum_dname,
            sampling_type="local",
            units="erg/s",
        )

        def _photon_emissivity_field(field, data):
            ret = data.ds.arr(self.process_data("photon_field", data),
                              "photons/s")
            return ret * data[ftype, "density"] / data[ftype, "mass"]

        phot_emiss_name = (
            ftype, emiss_name[1].replace("emissivity", "photon_emissivity")
        )

        ds.add_field(
            phot_emiss_name,
            function=_photon_emissivity_field,
            display_name=emiss_dname,
            sampling_type="local",
            units="photons/cm**3/s",
        )

        xray_fields = [emiss_name, lum_name, phot_emiss_name]

        if redshift > 0.0 or dist is not None:

            if dist is None:
                if cosmology is None:
                    if hasattr(ds, "cosmology"):
                        cosmology = ds.cosmology
                    else:
                        cosmology = Cosmology()
                D_L = cosmology.luminosity_distance(0.0, redshift)
                angular_scale = 1.0 / cosmology.angular_scale(0.0, redshift)
                dist_fac = ds.quan(
                    1.0 / (4.0 * np.pi * D_L * D_L * angular_scale * angular_scale).v,
                    "rad**-2",
                    )
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
                dist_fac = ds.quan(
                    1.0 / (4.0 * np.pi * dist * dist * angular_scale * angular_scale).v,
                    "rad**-2",
                    )

            ei_name = (
                ftype, emiss_name[1].replace("emissivity", "intensity")
            )
            ei_dname = emiss_dname.replace("\epsilon", "I")
            def _intensity_field(field, data):
                I = dist_fac * data[emiss_name]
                return I.in_units("erg/cm**3/s/arcsec**2")

            ds.add_field(
                ei_name,
                function=_intensity_field,
                display_name=ei_dname,
                sampling_type="local",
                units="erg/cm**3/s/arcsec**2",
            )

            i_name = (
                ftype, phot_emiss_name[1].replace("emissivity", "intensity")
            )
            i_dname = emiss_dname.replace("\epsilon", "I")

            def _photon_intensity_field(field, data):
                I = (1.0 + redshift) * dist_fac * data[phot_emiss_name]
                return I.in_units("photons/cm**3/s/arcsec**2")

            ds.add_field(
                i_name,
                function=_photon_intensity_field,
                display_name=i_dname,
                sampling_type="local",
                units="photons/cm**3/s/arcsec**2",
            )

            xray_fields += [ei_name, i_name]

        return xray_fields


metal_abund = {"angr": 0.0189,
               "aspl": 0.0134,
               "wilm": 0.0134,
               "lodd": 0.0133}


class ThermalSourceModel(SourceModel):
    def __init__(self, kT_min=0.025, kT_max=64.0, var_elem=None, 
                 max_density=5.0e-25, method="invert_cdf", 
                 abund_table="angr", prng=None):
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
        self.max_density = max_density
        self.tot_num_cells = 0  # Will be determined later
        self.ftype = "gas"
        self.abund_table = abund_table
        self.method = method
        self.prng = parse_prng(prng)
        self.kT_min = kT_min
        self.kT_max = kT_max

    def make_spectrum(self, data_source):
        self.setup_model(data_source, 0.0, 1.0)
        spec = np.zeros(self.spectral_model.emid.size)
        for chunk in data_source.chunks([], "io"):
            spec += self.process_data("spectrum", chunk)
        ebins = YTArray(self.spectral_model.ebins, "keV")
        return ebins, YTArray(spec, "photons/s")

    def process_data(self, mode, chunk):

        emid = self.spectral_model.emid
        ebins = self.spectral_model.ebins
        nchan = len(emid)

        orig_shape = chunk[self.temperature_field].shape
        if len(orig_shape) == 0:
            orig_ncells = 0
        else:
            orig_ncells = np.prod(orig_shape)
        if orig_ncells == 0:
            return

        ret = np.zeros(orig_ncells)

        cut = True

        if self.max_density is not None:
            cut &= np.ravel(chunk[self.density_field]) < self.max_density
        kT = np.ravel(
            chunk[self.temperature_field].to_value("keV", "thermal"))
        cut &= (kT >= self.kT_min) & (kT <= self.kT_max)

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

        if self.nei:
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
            elemZ = np.zeros((num_cells, self.num_var_elem))
            for j, key in enumerate(elem_keys):
                value = self.var_elem[key]
                if isinstance(value, float):
                    elemZ[:, j] = value
                else:
                    eZ = chunk[value]
                    fac = self.mconvert[key]
                    if str(eZ.units) != "Zsun":
                        fac /= X_H
                    elemZ[:, j] = np.ravel(eZ.d[cut]*fac)

        num_photons_max = 10000000
        number_of_photons = np.zeros(num_cells, dtype="int64")
        energies = np.zeros(num_photons_max)

        start_e = 0
        end_e = 0

        spec = np.zeros(nchan)
        idxs = np.where(cut)[0]

        for ck in chunked(range(num_cells), 100):

            ibegin = ck[0]
            iend = ck[-1]+1
            nck = iend-ibegin
            self.pbar.update(nck)

            cspec, mspec, vspec = self.spectral_model.get_spectrum(kT[ibegin:iend])

            tot_spec = cspec
            tot_spec += metalZ[ibegin:iend, np.newaxis] * mspec
            if vspec is not None:
                for j in range(self.num_var_elem):
                    tot_spec += elemZ[ibegin:iend, j, np.newaxis]*vspec[j, :, :]

            if mode in ["photons", "photon_field"]:

                cell_norm = tot_spec.sum(axis=-1)*cell_nrm[ibegin:iend]

                if mode == "photons":

                    cell_n = ensure_numpy_array(self.prng.poisson(lam=cell_norm))

                    number_of_photons[ibegin:iend] = cell_n
                    end_e += int(cell_n.sum())

                    norm_factor = 1.0 / tot_spec.sum(axis=-1)
                    p = norm_factor[:,np.newaxis]*tot_spec
                    cp = np.insert(np.cumsum(p, axis=-1), 0, 0.0, axis=1)
                    ei = start_e
                    for icell in range(nck):
                        cn = cell_n[icell]
                        if cn == 0:
                            continue
                        if self.method == "invert_cdf":
                            randvec = self.prng.uniform(size=cn)
                            randvec.sort()
                            cell_e = np.interp(randvec, cp[icell,:], ebins)
                        elif self.method == "accept_reject":
                            eidxs = self.prng.choice(nchan, size=cn, p=p[icell,:])
                            cell_e = emid[eidxs]
                        while ei+cn > num_photons_max:
                            num_photons_max *= 2
                        if num_photons_max > energies.size:
                            energies.resize(num_photons_max, refcheck=False)
                        energies[ei:ei+cn] = cell_e
                        ei += cn

                    start_e = end_e

                elif mode == "photon_field":

                    ret[idxs[ibegin:iend]] = cell_norm

            elif mode == "energy_field":

                ret[idxs[ibegin:iend]] = np.sum(
                    tot_spec*emid, axis=-1)*cell_nrm[ibegin:iend]

            elif mode == "spectrum":

                spec += np.sum(tot_spec*cell_nrm[ibegin:iend,np.newaxis], axis=0)

        if mode == "photons":
            active_cells = number_of_photons > 0
            idxs = idxs[active_cells]
            ncells = idxs.size
            return ncells, number_of_photons[active_cells], idxs, energies[:end_e].copy()
        elif mode == "spectrum":
            return spec
        else:
            return np.resize(ret, orig_shape)


class AtableSourceModel(ThermalSourceModel):
    def __init__(self, filename, temperature_field=None, norm_field=None, 
                 h_fraction=None, kT_min=0.025, kT_max=64.0, max_density=5.0e-25, 
                 var_elem=None, method="invert_cdf", abund_table="angr", prng=None):
        super().__init__(kT_min=kT_min, kT_max=kT_max, var_elem=var_elem, 
                         max_density=max_density, method=method, abund_table=abund_table, 
                         prng=prng)
        self.spectral_model = XSpecAtableModel(filename)
        self.temperature_field = temperature_field
        self.norm_field = norm_field
        self.h_fraction = h_fraction

class CIESourceModel(ThermalSourceModel):
    nei = False
    r"""
    Initialize a source model from a thermal spectrum.

    Parameters
    ----------
    emin : float
        The minimum energy for the spectrum in keV.
    emax : float
        The maximum energy for the spectrum in keV.
    nchan : integer
        The number of channels in the spectrum.
    Zmet : float, string, or tuple of strings
        The metallicity. If a float, assumes a constant metallicity throughout
        in solar units. If a string or tuple of strings, is taken to be the 
        name of the metallicity field.
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
        is taken to be the name of the hydrogen fraction field. Default: 0.74
    kT_min : float, optional
        The default minimum temperature in keV to compute emission for.
        Default: 0.025
    kT_max : float, optional
        The default maximum temperature in keV to compute emission for.
        Default: 64.0
    max_density : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
        The maximum density of the cells or particles to use when generating 
        photons. If a float, the units are assumed to be g/cm**3. 
        Default: 5e-25 g/cm**3.
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
    thermal_broad : boolean, optional
        Whether or not the spectral lines should be thermally
        broadened. Default: True
    model_root : string, optional
        The directory root where the model files are stored. If not provided,
        a default location known to pyXSIM is used. 
    model_vers : string, optional
        The version identifier string for the model files, e.g.
        "2.0.2". Default depends on the model used.
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
    prng : integer or :class:`~numpy.random.RandomState` object 
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, 
        such as for a test. Default is to use the :mod:`numpy.random` module.

    Examples
    --------
    >>> source_model = ThermalSourceModel("apec", 0.1, 10.0, 10000, 
    ...                                   ("gas", "metallicity"))
    """
    def __init__(self, emin, emax, nchan, Zmet, temperature_field=None, 
                 emission_measure_field=None, h_fraction=None, kT_min=0.025, 
                 kT_max=64.0, max_density=5.0e-25, var_elem=None, method="invert_cdf", 
                 thermal_broad=True, model_root=None, model_vers=None, nolines=False, 
                 abund_table="angr", prng=None):
        super().__init__(kT_min=kT_min, kT_max=kT_max, var_elem=var_elem,
                         max_density=max_density, method=method, abund_table=abund_table, 
                         prng=prng)
        self.temperature_field = temperature_field
        self.Zmet = Zmet
        self.emin = parse_value(emin, "keV")
        self.emax = parse_value(emax, "keV")
        self.spectral_model = TableApecModel(emin, emax, nchan,
                                             var_elem=self.var_elem_keys,
                                             thermal_broad=thermal_broad,
                                             model_root=model_root,
                                             model_vers=model_vers,
                                             nolines=nolines, nei=self.nei,
                                             abund_table=abund_table)
        self.var_elem_keys = self.spectral_model.var_elem_names
        self.var_ion_keys = self.spectral_model.var_ion_names
        self.spectral_norm = None
        self.redshift = None
        self.pbar = None
        self.emission_measure_field = emission_measure_field
        self.Zconvert = 1.0
        self.atable = self.spectral_model.atable
        self.mconvert = {}
        self.density_field = None  # Will be determined later
        if h_fraction is None:
            h_fraction = solar_H_abund
        self.h_fraction = h_fraction

    def setup_model(self, data_source, redshift, spectral_norm):
        if isinstance(data_source, Dataset):
            ds = data_source
        else:
            ds = data_source.ds
        if self.emission_measure_field is None:
            self.emission_measure_field = \
                (self.ftype, 'emission_measure')
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
        if not self.nei and not isinstance(self.Zmet, float):
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
        if self.temperature_field is None:
            self.temperature_field = (ftype, 'temperature')
        mylog.info(f"Using temperature field "
                   f"'{self.temperature_field}'.")
        self.spectral_model.prepare_spectrum(redshift, self.kT_min,
                                             self.kT_max)
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

    def cleanup_model(self):
        self.emission_measure_field = None
        self.temperature_field = None
        self.pbar.close()

    """
    def process_data(self, mode, chunk):

        emid = self.spectral_model.emid
        ebins = self.spectral_model.ebins
        nchan = len(emid)

        orig_shape = chunk[self.temperature_field].shape
        if len(orig_shape) == 0:
            orig_ncells = 0
        else:
            orig_ncells = np.prod(orig_shape)
        if orig_ncells == 0:
            return

        ret = np.zeros(orig_ncells)

        if self.max_density is None:
            dens_cut = 1
        else:
            dens_cut = chunk[self.density_field] < self.max_density
        kT = np.atleast_1d(
            chunk[self.temperature_field].to_value("keV", "thermal")).flat
        norm = np.atleast_1d(chunk[self.emission_measure_field].d*dens_cut).flat

        idxs = np.argsort(kT)

        kT_sorted = kT[idxs]
        idx_min = np.searchsorted(kT_sorted, self.kT_min)
        idx_max = np.searchsorted(kT_sorted, self.kT_max)
        idxs = idxs[idx_min:idx_max]
        num_cells = len(idxs)

        if mode == "photons":
            if num_cells == 0:
                self.pbar.update(orig_ncells)
                return
            else:
                self.pbar.update(orig_ncells-num_cells)
        elif num_cells == 0:
            return np.zeros(orig_shape)

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

        cell_nrm = norm[idxs]*self.spectral_norm

        if isinstance(self.h_fraction, float):
            X_H = self.h_fraction
        else:
            X_H = chunk[self.h_fraction].d[idxs]

        if self.nei:
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
                metalZ = np.atleast_1d(mZ.d[idxs]*fac).flat

        elemZ = None
        if self.num_var_elem > 0:
            elemZ = np.zeros((self.num_var_elem, num_cells))
            for j, key in enumerate(elem_keys):
                value = self.var_elem[key]
                if isinstance(value, float):
                    elemZ[j, :] = value
                else:
                    eZ = chunk[value]
                    fac = self.mconvert[key]
                    if str(eZ.units) != "Zsun":
                        fac /= X_H
                    elemZ[j, :] = np.atleast_1d(eZ.d[idxs]*fac).flat

        num_photons_max = 10000000
        number_of_photons = np.zeros(num_cells, dtype="int64")
        energies = np.zeros(num_photons_max)

        start_e = 0
        end_e = 0

        spec = np.zeros(nchan)

        for ibegin, iend, bcount, ikT in zip(bcell, ecell, bcounts, kT_idxs):

            self.pbar.update(bcount)

            kT = self.kT_bins[ikT] + 0.5*self.dkT[ikT]

            cnm = cell_nrm[ibegin:iend]

            cspec, mspec, vspec = self.spectral_model.get_spectrum(kT)

            if mode in ["photons", "photon_field"]:

                tot_ph_c = cspec.d.sum()
                tot_ph_m = mspec.d.sum()

                cell_norm_c = tot_ph_c*cnm
                cell_norm_m = tot_ph_m*metalZ[ibegin:iend]*cnm
                cell_norm = cell_norm_c + cell_norm_m

                if vspec is not None:
                    cell_norm_v = np.zeros(cnm.size)
                    for j in range(self.num_var_elem):
                        tot_ph_v = vspec.d[j, :].sum()
                        cell_norm_v += tot_ph_v*elemZ[j, ibegin:iend]*cnm
                    cell_norm += cell_norm_v

                if mode == "photons":

                    cell_n = ensure_numpy_array(self.prng.poisson(lam=cell_norm))

                    number_of_photons[ibegin:iend] = cell_n

                    end_e += int(cell_n.sum())

                    if self.method == "invert_cdf":
                        cumspec_c = np.pad(np.cumsum(cspec.d), (1,0))
                        cumspec_m = np.pad(np.cumsum(mspec.d), (1,0))
                        if vspec is None:
                            cumspec_v = None
                        else:
                            cumspec_v = np.zeros((self.num_var_elem, nchan+1))
                            for j in range(self.num_var_elem):
                                cumspec_v[j, 1:] = np.cumsum(vspec.d[j, :])

                    ei = start_e
                    for icell in range(ibegin, iend):
                        cn = number_of_photons[icell]
                        if cn == 0:
                            continue
                        # The rather verbose form of the few next statements is a
                        # result of code optimization and shouldn't be changed
                        # without checking for perfomance degradation.
                        if self.method == "invert_cdf":
                            cumspec = cumspec_c
                            cumspec += metalZ[icell] * cumspec_m
                            if cumspec_v is not None:
                                for j in range(self.num_var_elem):
                                    cumspec += elemZ[j, icell]*cumspec_v[j, :]
                            norm_factor = 1.0 / cumspec[-1]
                            cumspec *= norm_factor
                            randvec = self.prng.uniform(size=cn)
                            randvec.sort()
                            cell_e = np.interp(randvec, cumspec, ebins)
                        elif self.method == "accept_reject":
                            tot_spec = cspec.d
                            tot_spec += metalZ[icell] * mspec.d
                            if vspec is not None:
                                for j in range(self.num_var_elem):
                                    tot_spec += elemZ[j, icell]*vspec.d[j, :]
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

                elif mode == "photon_field":

                    ret[idxs[ibegin:iend]] = cell_norm

            elif mode == "energy_field":

                tot_e_c = (cspec.d*emid).sum()
                tot_e_m = (mspec.d*emid).sum()

                cell_e_c = tot_e_c * cnm
                cell_e_m = tot_e_m * metalZ[ibegin:iend] * cnm
                cell_e = cell_e_c + cell_e_m

                if vspec is not None:
                    cell_e_v = np.zeros(cnm.size)
                    for j in range(self.num_var_elem):
                        tot_e_v = (vspec.d[j, :]*emid).sum()
                        cell_e_v += tot_e_v * elemZ[j, ibegin:iend] * cnm
                    cell_e += cell_e_v

                ret[idxs[ibegin:iend]] = cell_e

            elif mode == "spectrum":

                spec += cnm.sum()*cspec.d
                spec += np.sum(metalZ[ibegin:iend]*cnm)*mspec.d
                if vspec is not None:
                    for j in range(self.num_var_elem):
                        spec += np.sum(elemZ[j, ibegin:iend]*cnm)*vspec.d[j, :]

        if mode == "photons":
            active_cells = number_of_photons > 0
            idxs = idxs[active_cells]
            ncells = idxs.size
            return ncells, number_of_photons[active_cells], idxs, energies[:end_e].copy()
        elif mode == "spectrum":
            return spec
        else:
            return np.resize(ret, orig_shape)
    """


class NEISourceModel(CIESourceModel):
    nei = True
    def __init__(self, emin, emax, nchan, var_elem, temperature_field=None,
                 emission_measure_field=None, h_fraction=None, kT_min=0.025,
                 kT_max=64.0, max_density=5.0e-25, method="invert_cdf", 
                 thermal_broad=True, model_root=None, model_vers=None, nolines=False, 
                 abund_table="angr", prng=None):
        super().__init__(emin, emax, nchan, 0.0, temperature_field=temperature_field,
                         emission_measure_field=emission_measure_field, h_fraction=h_fraction,
                         kT_min=kT_min, kT_max=kT_max, max_density=max_density, 
                         var_elem=var_elem, method=method, thermal_broad=thermal_broad, 
                         model_root=model_root, model_vers=model_vers, nolines=nolines, 
                         abund_table=abund_table, prng=prng)


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

    def make_spectrum(self, data_source, emin, emax, nbins):
        ebins = np.linspace(emin, emax, nbins+1)
        emid = 0.5*(ebins[1:]+ebins[:-1])
        spec = np.zeros(nbins)
        for chunk in data_source.chunks([], "io"):
            spec += self.process_data("spectrum", chunk, emid=emid)
        ebins = YTArray(ebins, "keV")
        return ebins, YTArray(spec, "photons/s")

    def process_data(self, mode, chunk, emid=None):

        num_cells = len(chunk[self.emission_field])

        if isinstance(self.alpha, float):
            alpha = self.alpha*np.ones(num_cells)
        else:
            alpha = chunk[self.alpha].v

        if mode in ["photons", "photon_field"]:

            norm_fac = self.emax.v**(1.-alpha) - self.emin.v**(1.-alpha)
            norm_fac[alpha == 1] = np.log(self.emax.v / self.emin.v)
            norm_fac *= self.e0.v**alpha
            norm = norm_fac * chunk[self.emission_field].v
            norm[alpha != 1] /= (1.-alpha[alpha != 1])

            if mode == "photons":

                norm *= self.spectral_norm*self.scale_factor

                number_of_photons = self.prng.poisson(lam=norm)

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
                        energies[start_e:end_e] = e * self.scale_factor
                        start_e = end_e

                active_cells = number_of_photons > 0
                ncells = active_cells.sum()

                return ncells, number_of_photons[active_cells], active_cells, energies[:end_e].copy()

            elif mode == "photon_field":

                return norm

        elif mode == "energy_field":

            norm_fac = self.emax.v**(2.-alpha) - self.emin.v**(2.-alpha)
            norm_fac *= self.e0.v ** alpha / (2. - alpha)

            return norm_fac * chunk[self.emission_field].v

        elif mode == "spectrum":

            return chunk[self.emission_field].v.sum()*(emid/self.e0.v)**(-alpha)


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
        If set to None (the default), it is assumed that the line is unbroadened.
    prng : integer or :class:`~numpy.random.RandomState` object 
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is to use the :mod:`numpy.random` module.

    Examples
    --------
    >>> location = (3.5, "keV")
    >>> sigma = (1000., "km/s")
    >>> line_model = LineEmissionSourceModel(location, "dark_matter_density_squared", sigma=sigma)
    """
    def __init__(self, e0, emission_field, sigma=None, prng=None):
        from unyt.exceptions import UnitConversionError
        self.e0 = parse_value(e0, "keV")
        if isinstance(sigma, Number):
            self.sigma = parse_value(sigma, "keV")
        elif isunitful(sigma):
            # The broadening is constant
            try:
                self.sigma = parse_value(sigma, "km/s")
                self.sigma *= self.e0/clight
                self.sigma.convert_to_units("keV")
            except UnitConversionError:
                self.sigma = parse_value(sigma, "keV")
        else:
            # Either no broadening or a field name
            self.sigma = sigma
        self.emission_field = emission_field
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

    def make_spectrum(self, data_source, emin, emax, nbins):
        ebins = np.linspace(emin, emax, nbins+1)
        spec = np.zeros(nbins)
        for chunk in data_source.chunks([], "io"):
            spec += self.process_data("spectrum", chunk, ebins=ebins)
        ebins = YTArray(ebins, "keV")
        return ebins, YTArray(spec, "photons/s")

    def process_data(self, mode, chunk, ebins=None):

        num_cells = len(chunk[self.emission_field])

        if mode == "photons":

            F = chunk[self.emission_field]*self.spectral_norm*self.scale_factor
            number_of_photons = self.prng.poisson(lam=F.in_cgs().d)

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

            energies = energies * self.scale_factor

            active_cells = number_of_photons > 0
            ncells = active_cells.sum()

            return ncells, number_of_photons[active_cells], active_cells, energies

        elif mode == "photon_field":

            return chunk[self.emission_field]

        elif mode == "energy_field":

            return self.e0*chunk[self.emission_field]

        elif mode == "spectrum":

            de = ebins-self.e0.value

            if isinstance(self.sigma, YTQuantity):
                xtmp = de/self.sigma.value
                ret = np.interp(xtmp, gx, gcdf)
                spec = chunk[self.emission_field].d.sum()*(ret[1:]-ret[:-1])
            elif self.sigma is not None:
                spec = np.zeros(ebins.size-1)
                sigma = (chunk[self.sigma]*self.e0/clight).to_value("keV")
                for i in range(num_cells):
                    xtmp = de/sigma[i]
                    ret = np.interp(xtmp, gx, gcdf)
                    spec += chunk[self.emission_field].d[i]*(ret[1:]-ret[:-1])
            else:
                spec = np.zeros(ebins.size-1)
                idx = np.searchsorted(ebins, self.e0.value)
                spec[idx] = chunk[self.emission_field].d.sum()

            return spec
