"""
Classes for generating lists of photons
"""
from six import string_types
from collections import defaultdict
import numpy as np
from yt.funcs import iterable
from pyxsim.utils import mylog, pixel_to_cel
from yt.utilities.physical_constants import clight
from yt.utilities.cosmology import Cosmology
from yt.utilities.orientation import Orientation
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    communication_system, get_mpi_type, parallel_capable, parallel_objects
from yt.units.yt_array import YTQuantity, YTArray, uconcatenate
import h5py
from pyxsim.spectral_models import absorb_models
from pyxsim.utils import parse_value, force_unicode, validate_parameters, \
    key_warning, ParameterDict
from pyxsim.event_list import EventList
from soxs.utils import parse_prng

comm = communication_system.communicators[-1]

axes_lookup = {"x": ("y","z"),
               "y": ("z","x"),
               "z": ("x","y")}

photon_units = {"energy": "keV",
                "dx": "kpc"}
for ax in "xyz":
    photon_units[ax] = "kpc"
    photon_units["v"+ax] = "km/s"

old_photon_keys = {"Energy": "energy",
                   "NumberOfPhotons": "num_photons"}
old_parameter_keys = {"FiducialExposureTime": "fid_exp_time",
                      "FiducialArea": "fid_area",
                      "FiducialRedshift": "fid_redshift",
                      "FiducialAngularDiameterDistance": "fid_d_a",
                      "HubbleConstant": "hubble",
                      "OmegaLambda": "omega_lambda",
                      "OmegaMatter": "omega_matter",
                      "DataType": "data_type"}

def determine_fields(ds, source_type):
    ds_type = ds.index.__class__.__name__
    if "ParticleIndex" in ds_type:
        position_fields = [(source_type, "particle_position_%s" % ax) for ax in "xyz"]
        velocity_fields = [(source_type, "particle_velocity_%s" % ax) for ax in "xyz"]
        width_field = (source_type, "smoothing_length")
    else:
        position_fields = [("index", ax) for ax in "xyz"]
        velocity_fields = [(source_type, "velocity_%s" % ax) for ax in "xyz"]
        width_field = ("index", "dx")
    if width_field not in ds.field_info:
        width_field = None
    return position_fields, velocity_fields, width_field

def concatenate_photons(photons):
    for key in photons:
        if len(photons[key]) > 0:
            photons[key] = uconcatenate(photons[key])
        elif key == "num_photons":
            photons[key] = np.array([])
        else:
            photons[key] = YTArray([], photon_units[key])

def find_object_bounds(data_source):
    # This logic is required to determine the bounds of 
    # the object, which is solely for fixing coordinates 
    # at periodic boundaries

    if hasattr(data_source, "base_object"):
        # This a cut region so we'll figure out
        # its bounds from its parent object
        data_src = data_source.base_object
    else:
        data_src = data_source

    if hasattr(data_src, "left_edge"):
        # Region or grid
        c = 0.5 * (data_src.left_edge + data_src.right_edge)
        w = data_src.right_edge - data_src.left_edge
        le = -0.5 * w + c
        re = 0.5 * w + c
    elif hasattr(data_src, "radius") and not hasattr(data_src, "height"):
        # Sphere
        le = -data_src.radius + data_src.center
        re = data_src.radius + data_src.center
    else:
        # Not sure what to do with any other object yet, so just
        # return the domain edges and punt.
        mylog.warning("You are using a region that is not currently "
                      "supported for straddling periodic boundaries. "
                      "Check to make sure that your region does not "
                      "do so.")
        le = data_source.ds.domain_left_edge
        re = data_source.ds.domain_right_edge

    return le.to("kpc"), re.to("kpc")

class PhotonList(object):

    def __init__(self, photons, parameters, cosmo):
        self.photons = photons
        self.parameters = ParameterDict(parameters, "PhotonList", old_parameter_keys)
        self.cosmo = cosmo
        self.num_cells = len(photons["x"])

        p_bins = np.cumsum(photons["num_photons"])
        self.p_bins = np.insert(p_bins, 0, [np.int64(0)])

    def keys(self):
        return self.photons.keys()

    def items(self):
        ret = []
        for k, v in self.photons.items():
            if k == "energy":
                ret.append((k, self[k]))
            else:
                ret.append((k,v))
        return ret

    def values(self):
        ret = []
        for k, v in self.photons.items():
            if k == "energy":
                ret.append(self[k])
            else:
                ret.append(v)
        return ret

    def __getitem__(self, key):
        if key in old_photon_keys:
            k = old_photon_keys[key]
            mylog.warning(key_warning % ("PhotonList", k))
        else:
            k = key
        if k == "energy":
            return [self.photons["energy"][self.p_bins[i]:self.p_bins[i+1]]
                    for i in range(self.num_cells)]
        else:
            return self.photons[k]

    def __contains__(self, key):
        if key in old_photon_keys:
            mylog.warning(key_warning % ("PhotonList", old_photon_keys[key]))
            return True
        return key in self.photons

    def __iter__(self):
        return iter(self.photons)

    def __repr__(self):
        return self.photons.__repr__()

    def __add__(self, other):
        validate_parameters(self.parameters, other.parameters)
        for param in ["hubble_constant", "omega_matter", "omega_lambda",
                      "omega_curvature"]:
            v1 = getattr(self.cosmo, param)
            v2 = getattr(other.cosmo, param)
            check_equal = np.allclose(np.array(v1), np.array(v2), rtol=0.0, atol=1.0e-10)
            if not check_equal:
                raise RuntimeError("The values for the parameter '%s' in the two" % param +
                                   " cosmologies are not identical (%s vs. %s)!" % (v1, v2))
        photons = {}
        for item1, item2 in zip(self.photons.items(), other.photons.items()):
            k1, v1 = item1
            k2, v2 = item2
            photons[k1] = uconcatenate([v1,v2])
        return PhotonList(photons, self.parameters, self.cosmo)

    @classmethod
    def from_file(cls, filename):
        r"""
        Initialize a :class:`~pyxsim.photon_list.PhotonList` from the HDF5 file *filename*.
        """

        photons = {}
        parameters = {}

        f = h5py.File(filename, "r")

        p = f["/parameters"]
        parameters["fid_exp_time"] = YTQuantity(p["fid_exp_time"].value, "s")
        parameters["fid_area"] = YTQuantity(p["fid_area"].value, "cm**2")
        parameters["fid_redshift"] = p["fid_redshift"].value
        parameters["fid_d_a"] = YTQuantity(p["fid_d_a"].value, "Mpc")
        parameters["hubble"] = p["hubble"].value
        parameters["omega_matter"] = p["omega_matter"].value
        parameters["omega_lambda"] = p["omega_lambda"].value
        if "data_type" in p:
            parameters["data_type"] = force_unicode(p["data_type"].value)
        else:
            parameters["data_type"] = "cells"

        d = f["/data"]

        num_cells = d["x"].size
        start_c = comm.rank*num_cells//comm.size
        end_c = (comm.rank+1)*num_cells//comm.size

        photons["x"] = YTArray(d["x"][start_c:end_c], "kpc")
        photons["y"] = YTArray(d["y"][start_c:end_c], "kpc")
        photons["z"] = YTArray(d["z"][start_c:end_c], "kpc")
        photons["dx"] = YTArray(d["dx"][start_c:end_c], "kpc")
        photons["vx"] = YTArray(d["vx"][start_c:end_c], "km/s")
        photons["vy"] = YTArray(d["vy"][start_c:end_c], "km/s")
        photons["vz"] = YTArray(d["vz"][start_c:end_c], "km/s")

        n_ph = d["num_photons"][:]

        if comm.rank == 0:
            start_e = np.int64(0)
        else:
            start_e = n_ph[:start_c].sum()
        end_e = start_e + np.int64(n_ph[start_c:end_c].sum())

        photons["num_photons"] = n_ph[start_c:end_c]
        photons["energy"] = YTArray(d["energy"][start_e:end_e], "keV")

        f.close()

        cosmo = Cosmology(hubble_constant=parameters["hubble"],
                          omega_matter=parameters["omega_matter"],
                          omega_lambda=parameters["omega_lambda"])

        return cls(photons, parameters, cosmo)

    @classmethod
    def from_data_source(cls, data_source, redshift, area,
                         exp_time, source_model, parameters=None,
                         center=None, dist=None, cosmology=None,
                         velocity_fields=None):
        r"""
        Initialize a :class:`~pyxsim.photon_list.PhotonList` from a yt data source.
        The redshift, collecting area, exposure time, and cosmology are stored in the
        *parameters* dictionary which is passed to the *source_model* function.

        Parameters
        ----------
        data_source : :class:`~yt.data_objects.data_containers.YTSelectionContainer`
            The data source from which the photons will be generated.
        redshift : float
            The cosmological redshift for the photons.
        area : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The collecting area to determine the number of photons. If units are
            not specified, it is assumed to be in cm^2.
        exp_time : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The exposure time to determine the number of photons. If units are
            not specified, it is assumed to be in seconds.
        source_model : :class:`~pyxsim.source_models.SourceModel`
            A source model used to generate the photons.
        parameters : dict, optional
            A dictionary of parameters to be passed for the source model to use, if necessary.
        center : string or array_like, optional
            The origin of the photon spatial coordinates. Accepts "c", "max", or a coordinate. 
            If not specified, pyxsim attempts to use the "center" field parameter of the data_source. 
        dist : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The angular diameter distance, used for nearby sources. This may be
            optionally supplied instead of it being determined from the *redshift*
            and given *cosmology*. If units are not specified, it is assumed to be
            in Mpc. To use this, the redshift must be set to zero. 
        cosmology : :class:`~yt.utilities.cosmology.Cosmology`, optional
            Cosmological information. If not supplied, we try to get
            the cosmology from the dataset. Otherwise, LCDM with
            the default yt parameters is assumed.
        velocity_fields : list of fields
            The yt fields to use for the velocity. If not specified, the following will
            be assumed:
            ['velocity_x', 'velocity_y', 'velocity_z'] for grid datasets
            ['particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z'] for particle datasets

        Examples
        --------
        >>> thermal_model = ThermalSourceModel(apec_model, Zmet=0.3)
        >>> redshift = 0.05
        >>> area = 6000.0 # assumed here in cm**2
        >>> time = 2.0e5 # assumed here in seconds
        >>> sp = ds.sphere("c", (500., "kpc"))
        >>> my_photons = PhotonList.from_data_source(sp, redshift, area,
        ...                                          time, thermal_model)
        """
        ds = data_source.ds

        if parameters is None:
             parameters = {}
        if cosmology is None:
            if hasattr(ds, 'cosmology'):
                cosmo = ds.cosmology
            else:
                cosmo = Cosmology()
        else:
            cosmo = cosmology
        mylog.info("Cosmology: h = %g, omega_matter = %g, omega_lambda = %g" %
                   (cosmo.hubble_constant, cosmo.omega_matter, cosmo.omega_lambda))
        if dist is None:
            if redshift <= 0.0:
                msg = "If redshift <= 0.0, you must specify a distance to the source using the 'dist' argument!"
                mylog.error(msg)
                raise ValueError(msg)
            D_A = cosmo.angular_diameter_distance(0.0, redshift).in_units("Mpc")
        else:
            D_A = parse_value(dist, "Mpc")
            if redshift > 0.0:
                mylog.warning("Redshift must be zero for nearby sources. Resetting redshift to 0.0.")
                redshift = 0.0

        if isinstance(center, string_types):
            if center == "center" or center == "c":
                parameters["center"] = ds.domain_center
            elif center == "max" or center == "m":
                parameters["center"] = ds.find_max("density")[-1]
        elif iterable(center):
            if isinstance(center, YTArray):
                parameters["center"] = center.in_units("code_length")
            elif isinstance(center, tuple):
                if center[0] == "min":
                    parameters["center"] = ds.find_min(center[1])[-1]
                elif center[0] == "max":
                    parameters["center"] = ds.find_max(center[1])[-1]
                else:
                    raise RuntimeError
            else:
                parameters["center"] = ds.arr(center, "code_length")
        elif center is None:
            if hasattr(data_source, "left_edge"):
                parameters["center"] = 0.5*(data_source.left_edge+data_source.right_edge)
            else:
                parameters["center"] = data_source.get_field_parameter("center")

        parameters["fid_exp_time"] = parse_value(exp_time, "s")
        parameters["fid_area"] = parse_value(area, "cm**2")
        parameters["fid_redshift"] = redshift
        parameters["fid_d_a"] = D_A
        parameters["hubble"] = cosmo.hubble_constant
        parameters["omega_matter"] = cosmo.omega_matter
        parameters["omega_lambda"] = cosmo.omega_lambda

        D_A = parameters["fid_d_a"].in_cgs()
        dist_fac = 1.0/(4.*np.pi*D_A.value*D_A.value*(1.+redshift)**2)
        spectral_norm = parameters["fid_area"].v*parameters["fid_exp_time"].v*dist_fac

        source_model.setup_model(data_source, redshift, spectral_norm)

        p_fields, v_fields, w_field = determine_fields(ds, source_model.source_type)

        if velocity_fields is not None:
            v_fields = velocity_fields

        if p_fields[0] == ("index", "x"):
            parameters["data_type"] = "cells"
        else:
            parameters["data_type"] = "particles"

        citer = data_source.chunks([], "io")

        photons = defaultdict(list)

        for chunk in parallel_objects(citer):

            chunk_data = source_model(chunk)

            if chunk_data is not None:
                number_of_photons, idxs, energies = chunk_data
                photons["num_photons"].append(number_of_photons)
                photons["energy"].append(ds.arr(energies, "keV"))
                photons["x"].append(chunk[p_fields[0]][idxs].in_units("kpc"))
                photons["y"].append(chunk[p_fields[1]][idxs].in_units("kpc"))
                photons["z"].append(chunk[p_fields[2]][idxs].in_units("kpc"))
                photons["vx"].append(chunk[v_fields[0]][idxs].in_units("km/s"))
                photons["vy"].append(chunk[v_fields[1]][idxs].in_units("km/s"))
                photons["vz"].append(chunk[v_fields[2]][idxs].in_units("km/s"))
                if w_field is None:
                    photons["dx"].append(ds.arr(np.zeros(len(idxs)), "kpc"))
                else:
                    photons["dx"].append(chunk[w_field][idxs].in_units("kpc"))

        source_model.cleanup_model()

        concatenate_photons(photons)

        c = parameters["center"].to("kpc")

        if sum(ds.periodicity) > 0:
            # Fix photon coordinates for regions crossing a periodic boundary
            dw = ds.domain_width.to("kpc")
            le, re = find_object_bounds(data_source)
            for i, ax in enumerate("xyz"):
                if ds.periodicity[i] and len(photons[ax]) > 0:
                    tfl = photons[ax] < le[i]
                    tfr = photons[ax] > re[i]
                    photons[ax][tfl] += dw[i]
                    photons[ax][tfr] -= dw[i]

        # Re-center all coordinates
        for i, ax in enumerate("xyz"):
            if len(photons[ax]) > 0:
                photons[ax] -= c[i]

        mylog.info("Finished generating photons.")
        mylog.info("Number of photons generated: %d" % int(np.sum(photons["num_photons"])))
        mylog.info("Number of cells with photons: %d" % len(photons["x"]))

        return cls(photons, parameters, cosmo)

    def write_h5_file(self, photonfile):
        """
        Write the :class:`~pyxsim.photon_list.PhotonList` to the HDF5 file *photonfile*.
        """

        if parallel_capable:

            mpi_long = get_mpi_type("int64")
            mpi_double = get_mpi_type("float64")

            local_num_cells = len(self.photons["x"])
            sizes_c = comm.comm.gather(local_num_cells, root=0)

            local_num_photons = np.sum(self.photons["num_photons"])
            sizes_p = comm.comm.gather(local_num_photons, root=0)

            if comm.rank == 0:
                num_cells = sum(sizes_c)
                num_photons = sum(sizes_p)
                disps_c = [sum(sizes_c[:i]) for i in range(len(sizes_c))]
                disps_p = [sum(sizes_p[:i]) for i in range(len(sizes_p))]
                x = np.zeros(num_cells)
                y = np.zeros(num_cells)
                z = np.zeros(num_cells)
                vx = np.zeros(num_cells)
                vy = np.zeros(num_cells)
                vz = np.zeros(num_cells)
                dx = np.zeros(num_cells)
                n_ph = np.zeros(num_cells, dtype="int64")
                e = np.zeros(num_photons)
            else:
                sizes_c = []
                sizes_p = []
                disps_c = []
                disps_p = []
                x = np.empty([])
                y = np.empty([])
                z = np.empty([])
                vx = np.empty([])
                vy = np.empty([])
                vz = np.empty([])
                dx = np.empty([])
                n_ph = np.empty([])
                e = np.empty([])

            comm.comm.Gatherv([self.photons["x"].d, local_num_cells, mpi_double],
                              [x, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["y"].d, local_num_cells, mpi_double],
                              [y, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["z"].d, local_num_cells, mpi_double],
                              [z, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["vx"].d, local_num_cells, mpi_double],
                              [vx, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["vy"].d, local_num_cells, mpi_double],
                              [vy, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["vz"].d, local_num_cells, mpi_double],
                              [vz, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["dx"].d, local_num_cells, mpi_double],
                              [dx, (sizes_c, disps_c), mpi_double], root=0)
            comm.comm.Gatherv([self.photons["num_photons"], local_num_cells, mpi_long],
                              [n_ph, (sizes_c, disps_c), mpi_long], root=0)
            comm.comm.Gatherv([self.photons["energy"].d, local_num_photons, mpi_double],
                              [e, (sizes_p, disps_p), mpi_double], root=0)

        else:

            x = self.photons["x"].d
            y = self.photons["y"].d
            z = self.photons["z"].d
            vx = self.photons["vx"].d
            vy = self.photons["vy"].d
            vz = self.photons["vz"].d
            dx = self.photons["dx"].d
            n_ph = self.photons["num_photons"]
            e = self.photons["energy"].d

        if comm.rank == 0:

            f = h5py.File(photonfile, "w")

            # Parameters

            p = f.create_group("parameters")
            p.create_dataset("fid_area", data=float(self.parameters["fid_area"]))
            p.create_dataset("fid_exp_time", data=float(self.parameters["fid_exp_time"]))
            p.create_dataset("fid_redshift", data=self.parameters["fid_redshift"])
            p.create_dataset("hubble", data=self.parameters["hubble"])
            p.create_dataset("omega_matter", data=self.parameters["omega_matter"])
            p.create_dataset("omega_lambda", data=self.parameters["omega_lambda"])
            p.create_dataset("fid_d_a", data=float(self.parameters["fid_d_a"]))
            p.create_dataset("data_type", data=self.parameters["data_type"])

            # Data

            d = f.create_group("data")
            d.create_dataset("x", data=x)
            d.create_dataset("y", data=y)
            d.create_dataset("z", data=z)
            d.create_dataset("vx", data=vx)
            d.create_dataset("vy", data=vy)
            d.create_dataset("vz", data=vz)
            d.create_dataset("dx", data=dx)
            d.create_dataset("num_photons", data=n_ph)
            d.create_dataset("energy", data=e)

            f.close()

        comm.barrier()

    def project_photons(self, normal, sky_center, area_new=None, 
                        exp_time_new=None, redshift_new=None, 
                        dist_new=None, absorb_model=None, nH=None,
                        no_shifting=False, north_vector=None,
                        smooth_positions=None, prng=None):
        r"""
        Projects photons onto an image plane given a line of sight.
        Returns a new :class:`~pyxsim.event_list.EventList`.

        Parameters
        ----------
        normal : character or array-like
            Normal vector to the plane of projection. If "x", "y", or "z", will
            assume to be along that axis (and will probably be faster). Otherwise,
            should be an off-axis normal vector, e.g [1.0, 2.0, -3.0]
        sky_center : array-like
            Center RA, Dec of the events in degrees.
        area_new : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            A value for the (constant) collecting area of the detector. If
            units are not specified, is assumed to be in cm**2.
        exp_time_new : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            A new value for the exposure time. If units are not specified
            it is assumed to be in seconds.
        redshift_new : float, optional
            A new value for the cosmological redshift. Cannot be specified
            if you applied foreground galactic absorption already in the 
            :class:`~pyxsim.photon_list.PhotonList` instance.
        dist_new : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The new value for the angular diameter distance, used for nearby sources.
            This may be optionally supplied instead of it being determined from the
            cosmology. If units are not specified, it is assumed to be in Mpc. To 
            use this, the redshift must be zero. 
        absorb_model : string or :class:`~pyxsim.spectral_models.AbsorptionModel` 
            A model for foreground galactic absorption, to simulate the absorption
            of events before being detected. This cannot be applied here if you 
            already did this step previously in the creation of the 
            :class:`~pyxsim.photon_list.PhotonList` instance. Known options for 
            strings are "wabs" and "tbabs".
        nH : float, optional
            The foreground column density in units of 10^22 cm^{-2}. Only used if
            absorption is applied.
        no_shifting : boolean, optional
            If set, the photon energies will not be Doppler shifted.
        north_vector : a sequence of floats
            A vector defining the "up" direction. This option sets the orientation of
            the plane of projection. If not set, an arbitrary grid-aligned north_vector
            is chosen. Ignored in the case where a particular axis (e.g., "x", "y", or
            "z") is explicitly specified.
        smooth_positions : float, optional
            Apply a gaussian smoothing operation to the sky positions of the events. 
            This may be useful when the binned events appear blocky due to their uniform
            distribution within simulation cells. However, this will move the events away
            from their originating position on the sky, and so may distort surface brightness
            profiles and/or spectra. Should probably only be used for visualization purposes.
            Supply a float here to smooth with a standard deviation with this fraction 
            of the cell or particle size. Default: None
        prng : integer or :class:`~numpy.random.RandomState` object 
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is to use the :mod:`numpy.random` module.

        Examples
        --------
        >>> L = np.array([0.1,-0.2,0.3])
        >>> events = my_photons.project_photons(L, [30., 45.], area_new=10000.,
        ...                                     redshift_new=0.05)
        """
        prng = parse_prng(prng)

        if smooth_positions is not None and self.parameters["data_type"] == "particles":
            raise RuntimeError("The 'smooth_positions' argument should not be used with "
                               "particle-based datasets!")

        change_redshift = redshift_new is not None
        change_dist = dist_new is not None

        if change_redshift and change_dist:
            raise RuntimeError("You may specify a new redshift or distance, "
                               "but not both!")

        if isinstance(absorb_model, string_types):
            if absorb_model not in absorb_models:
                raise KeyError("%s is not a known absorption model!" % absorb_model)
            absorb_model = absorb_models[absorb_model]
        if absorb_model is not None:
            if nH is None:
                raise RuntimeError("You specified an absorption model, but didn't "
                                   "specify a value for nH!")
            absorb_model = absorb_model(nH)

        sky_center = YTArray(sky_center, "degree")

        dx = self.photons["dx"].d

        if not isinstance(normal, string_types):
            L = np.array(normal)
            orient = Orientation(L, north_vector=north_vector)
            x_hat = orient.unit_vectors[0]
            y_hat = orient.unit_vectors[1]
            z_hat = orient.unit_vectors[2]

        n_ph = self.photons["num_photons"]
        n_ph_tot = n_ph.sum()

        parameters = {}

        zobs0 = self.parameters["fid_redshift"]
        D_A0 = self.parameters["fid_d_a"]

        scale_factor = 1.0

        if (exp_time_new is None and area_new is None and
            redshift_new is None and dist_new is None):
            my_n_obs = n_ph_tot
            D_A = D_A0
        else:
            if exp_time_new is None:
                Tratio = 1.
            else:
                exp_time_new = parse_value(exp_time_new, "s")
                Tratio = exp_time_new/self.parameters["fid_exp_time"]
            if area_new is None:
                Aratio = 1.
            else:
                area_new = parse_value(area_new, "cm**2")
                Aratio = area_new/self.parameters["fid_area"]
            if redshift_new is None and dist_new is None:
                Dratio = 1.
                D_A = D_A0
            else:
                if dist_new is not None:
                    if redshift_new is not None and redshift_new > 0.0:
                        mylog.warning("Redshift must be zero for nearby sources. "
                                      "Assuming redshift of 0.0.")
                    D_A = parse_value(dist_new, "Mpc")
                else:
                    zobs = redshift_new
                    D_A = self.cosmo.angular_diameter_distance(0.0, zobs).in_units("Mpc")
                Dratio = D_A0*D_A0*(1.+zobs0)**3 / \
                         (D_A*D_A*(1.+zobs)**3)
                scale_factor = (1.+zobs0)/(1.+zobs)
            fak = Aratio*Tratio*Dratio
            if fak > 1:
                raise ValueError("This combination of requested parameters results in "
                                 "%g%% more photons collected than are " % (100.*(fak-1.)) +
                                 "available in the sample. Please reduce the collecting "
                                 "area, exposure time, or increase the distance/redshift "
                                 "of the object. Alternatively, generate a larger sample "
                                 "of photons.")
            my_n_obs = np.int64(n_ph_tot*fak)

        if my_n_obs == n_ph_tot:
            idxs = np.arange(my_n_obs, dtype='int64')
        else:
            idxs = prng.permutation(n_ph_tot)[:my_n_obs].astype("int64")
        obs_cells = np.searchsorted(self.p_bins, idxs, side='right')-1
        delta = dx[obs_cells]

        events = {}

        eobs = self.photons["energy"][idxs]

        if not no_shifting:
            if isinstance(normal, string_types):
                shift = -self.photons["v%s" % normal][obs_cells].in_cgs()
            else:
                shift = -(self.photons["vx"]*z_hat[0] +
                          self.photons["vy"]*z_hat[1] +
                          self.photons["vz"]*z_hat[2])[obs_cells].in_cgs()
            shift /= clight
            np.sqrt((1.-shift)/(1.+shift), shift)
            np.multiply(eobs, shift, eobs)

        eobs *= scale_factor

        if absorb_model is None:
            detected = np.ones(eobs.shape, dtype='bool')
        else:
            detected = absorb_model.absorb_photons(eobs, prng=prng)

        num_det = detected.sum()

        events["eobs"] = eobs[detected]

        if num_det > 0:

            deld = delta[detected]
            ocells = obs_cells[detected]

            if isinstance(normal, string_types):

                if self.parameters["data_type"] == "cells":
                    xsky = prng.uniform(low=-0.5, high=0.5, size=num_det)
                    ysky = prng.uniform(low=-0.5, high=0.5, size=num_det)
                elif self.parameters["data_type"] == "particles":
                    xsky = prng.normal(loc=0.0, scale=1.0, size=num_det)
                    ysky = prng.normal(loc=0.0, scale=1.0, size=num_det)

                np.multiply(xsky, deld, xsky)
                np.multiply(ysky, deld, ysky)
                np.add(xsky, self.photons[axes_lookup[normal][0]].d[ocells], xsky)
                np.add(ysky, self.photons[axes_lookup[normal][1]].d[ocells], ysky)

            else:

                if self.parameters["data_type"] == "cells":
                    r = prng.uniform(low=-0.5, high=0.5, size=(3, num_det))
                elif self.parameters["data_type"] == "particles":
                    r = prng.normal(loc=0.0, scale=1.0, size=(3, num_det))

                np.multiply(r, deld, r)
                r[0,:] += self.photons["x"].d[ocells]
                r[1,:] += self.photons["y"].d[ocells]
                r[2,:] += self.photons["z"].d[ocells]

                xsky, ysky = np.dot([x_hat, y_hat], r)

            if smooth_positions is not None:
                sigma = smooth_positions*deld
                xsky += sigma*prng.normal(loc=0.0, scale=1.0, size=num_det)
                ysky += sigma*prng.normal(loc=0.0, scale=1.0, size=num_det)

            d_a = D_A.to("kpc").v
            xsky /= d_a
            ysky /= d_a

            xsky, ysky = pixel_to_cel(xsky, ysky, sky_center)

        else:

            xsky = []
            ysky = []

        events["xsky"] = YTArray(xsky, "degree")
        events["ysky"] = YTArray(ysky, "degree")

        num_events = comm.mpi_allreduce(num_det)

        if comm.rank == 0:
            mylog.info("Total number of observed photons: %d" % num_events)

        if exp_time_new is None:
            parameters["exp_time"] = self.parameters["fid_exp_time"]
        else:
            parameters["exp_time"] = exp_time_new
        if area_new is None:
            parameters["area"] = self.parameters["fid_area"]
        else:
            parameters["area"] = area_new
        parameters["sky_center"] = sky_center

        return EventList(events, parameters)
