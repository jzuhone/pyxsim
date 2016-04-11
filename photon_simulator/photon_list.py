"""
Classes for generating lists of photons and detected events

The SciPy Proceeding that describes this module in detail may be found at:

http://conference.scipy.org/proceedings/scipy2014/zuhone.html

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
from collections import defaultdict
import numpy as np
from yt.funcs import mylog, get_pbar, iterable, ensure_list
from yt.utilities.physical_constants import clight
from yt.utilities.cosmology import Cosmology
from yt.utilities.orientation import Orientation
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    communication_system, get_mpi_type, parallel_capable, parallel_objects
from yt.units.yt_array import YTQuantity, YTArray, uconcatenate
from yt.utilities.on_demand_imports import _astropy
import h5py
from photon_simulator.utils import parse_value, force_unicode
from photon_simulator.event_list import EventList

comm = communication_system.communicators[-1]

axes_lookup = {"x": ("y","z"),
               "y": ("z","x"),
               "z": ("x","y")}

photon_units = {"Energy": "keV",
                "dx": "kpc"}
for ax in "xyz":
    photon_units[ax] = "kpc"
    photon_units["v"+ax] = "km/s"

def determine_fields(ds):
    ds_type = ds.index.__class__.__name__
    if "ParticleIndex" in ds_type:
        position_fields = ["particle_position_%s" % ax for ax in "xyz"]
        velocity_fields = ["particle_velocity_%s" % ax for ax in "xyz"]
        width_field = "smoothing_length"
    else:
        position_fields = list("xyz")
        velocity_fields = ["velocity_%s" % ax for ax in "xyz"]
        width_field = "dx"
    return position_fields, velocity_fields, width_field

def concatenate_photons(photons):
    for key in photons:
        if len(photons[key]) > 0:
            photons[key] = uconcatenate(photons[key])
        elif key == "NumberOfPhotons":
            photons[key] = np.array([])
        else:
            photons[key] = YTArray([], photon_units[key])

class PhotonList(object):

    def __init__(self, photons, parameters, cosmo):
        self.photons = photons
        self.parameters = parameters
        self.cosmo = cosmo
        self.num_cells = len(photons["x"])

        p_bins = np.cumsum(photons["NumberOfPhotons"])
        self.p_bins = np.insert(p_bins, 0, [np.uint64(0)])

    def keys(self):
        return self.photons.keys()

    def items(self):
        ret = []
        for k, v in self.photons.items():
            if k == "Energy":
                ret.append((k, self[k]))
            else:
                ret.append((k,v))
        return ret

    def values(self):
        ret = []
        for k, v in self.photons.items():
            if k == "Energy":
                ret.append(self[k])
            else:
                ret.append(v)
        return ret

    def __getitem__(self, key):
        if key == "Energy":
            return [self.photons["Energy"][self.p_bins[i]:self.p_bins[i+1]]
                    for i in range(self.num_cells)]
        else:
            return self.photons[key]

    def __contains__(self, key):
        return key in self.photons

    def __repr__(self):
        return self.photons.__repr__()

    @classmethod
    def from_file(cls, filename):
        r"""
        Initialize a PhotonList from the HDF5 file *filename*.
        """

        photons = {}
        parameters = {}

        f = h5py.File(filename, "r")

        p = f["/parameters"]
        parameters["FiducialExposureTime"] = YTQuantity(p["fid_exp_time"].value, "s")
        parameters["FiducialArea"] = YTQuantity(p["fid_area"].value, "cm**2")
        parameters["FiducialRedshift"] = p["fid_redshift"].value
        parameters["FiducialAngularDiameterDistance"] = YTQuantity(p["fid_d_a"].value, "Mpc")
        parameters["Dimension"] = p["dimension"].value
        parameters["Width"] = YTQuantity(p["width"].value, "kpc")
        parameters["HubbleConstant"] = p["hubble"].value
        parameters["OmegaMatter"] = p["omega_matter"].value
        parameters["OmegaLambda"] = p["omega_lambda"].value
        if "data_type" in p:
            parameters["DataType"] = force_unicode(p["data_type"].value)
        else:
            parameters["DataType"] = "cells"

        d = f["/data"]

        num_cells = d["x"][:].shape[0]
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
            start_e = np.uint64(0)
        else:
            start_e = n_ph[:start_c].sum()
        end_e = start_e + np.uint64(n_ph[start_c:end_c].sum())

        photons["NumberOfPhotons"] = n_ph[start_c:end_c]
        photons["Energy"] = YTArray(d["energy"][start_e:end_e], "keV")

        f.close()

        cosmo = Cosmology(hubble_constant=parameters["HubbleConstant"],
                          omega_matter=parameters["OmegaMatter"],
                          omega_lambda=parameters["OmegaLambda"])

        return cls(photons, parameters, cosmo)

    @classmethod
    def from_data_source(cls, data_source, redshift, area,
                         exp_time, source_model, parameters=None,
                         center=None, dist=None, cosmology=None,
                         velocity_fields=None):
        r"""
        Initialize a PhotonList from a photon model. The redshift, collecting area,
        exposure time, and cosmology are stored in the *parameters* dictionary which
        is passed to the *photon_model* function. 

        Parameters
        ----------
        data_source : `yt.data_objects.data_containers.YTSelectionContainer`
            The data source from which the photons will be generated. NOTE: The
            PointSourceModel does not require a data_source, so *None* can be
            supplied here in that case.
        redshift : float
            The cosmological redshift for the photons.
        area : float, (value, unit) tuple, or YTQuantity.
            The collecting area to determine the number of photons. If units are
            not specified, it is assumed to be in cm**2.
        exp_time : float, (value, unit) tuple, or YTQuantity.
            The exposure time to determine the number of photons. If units are
            not specified, it is assumed to be in seconds.
        source_model : function
            A function that takes the *data_source* and the *parameters*
            dictionary and returns a *photons* dictionary. Must be of the
            form: photon_model(data_source, parameters)
        parameters : dict, optional
            A dictionary of parameters to be passed to the user function.
        center : string or array_like, optional
            The origin of the photons. Accepts "c", "max", or a coordinate.
        dist : float, (value, unit) tuple, or YTQuantity, optional
            The angular diameter distance, used for nearby sources. This may be
            optionally supplied instead of it being determined from the *redshift*
            and given *cosmology*. If units are not specified, it is assumed to be
            in Mpc.
        cosmology : `yt.utilities.cosmology.Cosmology`, optional
            Cosmological information. If not supplied, we try to get
            the cosmology from the dataset. Otherwise, \LambdaCDM with
            the default yt parameters is assumed.

        Examples
        --------
        This is the simplest possible example, where we call the built-in thermal model:

        >>> thermal_model = ThermalSourceModel(apec_model, Zmet=0.3)
        >>> redshift = 0.05
        >>> area = 6000.0 # assumed here in cm**2
        >>> time = 2.0e5 # assumed here in seconds
        >>> sp = ds.sphere("c", (500., "kpc"))
        >>> my_photons = PhotonList.from_data_source(sp, redshift, area,
        ...                                          time, thermal_model)
        """
        ds = data_source.ds

        p_fields, v_fields, w_field = determine_fields(ds)

        if velocity_fields is not None:
            v_fields = velocity_fields

        if parameters is None:
             parameters = {}
        if cosmology is None:
            hubble = getattr(ds, "hubble_constant", None)
            omega_m = getattr(ds, "omega_matter", None)
            omega_l = getattr(ds, "omega_lambda", None)
            if hubble == 0:
                hubble = None
            if hubble is not None and \
               omega_m is not None and \
               omega_l is not None:
                cosmo = Cosmology(hubble_constant=hubble,
                                  omega_matter=omega_m,
                                  omega_lambda=omega_l)
            else:
                cosmo = Cosmology()
        else:
            cosmo = cosmology
        mylog.info("Cosmology: h = %g, omega_matter = %g, omega_lambda = %g" %
                   (cosmo.hubble_constant, cosmo.omega_matter, cosmo.omega_lambda))
        if dist is None:
            D_A = cosmo.angular_diameter_distance(0.0,redshift).in_units("Mpc")
        else:
            D_A = parse_value(dist, "Mpc")
            redshift = 0.0

        if center in ("center", "c"):
            parameters["center"] = ds.domain_center
        elif center in ("max", "m"):
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
            parameters["center"] = data_source.get_field_parameter("center")

        parameters["FiducialExposureTime"] = parse_value(exp_time, "s")
        parameters["FiducialArea"] = parse_value(area, "cm**2")
        parameters["FiducialRedshift"] = redshift
        parameters["FiducialAngularDiameterDistance"] = D_A
        parameters["HubbleConstant"] = cosmo.hubble_constant
        parameters["OmegaMatter"] = cosmo.omega_matter
        parameters["OmegaLambda"] = cosmo.omega_lambda

        if p_fields[0] == "x":
            parameters["DataType"] = "cells"
        else:
            parameters["DataType"] = "particles"

        dimension = 0
        width = 0.0
        for i, ax in enumerate("xyz"):
            le, re = data_source.quantities.extrema(ax)
            delta_min, delta_max = data_source.quantities.extrema("d%s"%ax)
            le -= 0.5*delta_max
            re += 0.5*delta_max
            width = max(width, re-parameters["center"][i], parameters["center"][i]-le)
            dimension = max(dimension, int(width/delta_min))
        parameters["Dimension"] = 2*dimension
        parameters["Width"] = 2.*width.in_units("kpc")

        D_A = parameters["FiducialAngularDiameterDistance"].in_cgs()
        dist_fac = 1.0/(4.*np.pi*D_A.value*D_A.value*(1.+redshift)**2)
        spectral_norm = parameters["FiducialArea"].v*parameters["FiducialExposureTime"].v*dist_fac

        citer = data_source.chunks([], "io")

        photons = defaultdict(list)

        source_model.setup_model(data_source, redshift, spectral_norm)

        for chunk in parallel_objects(citer):

            chunk_data = source_model(chunk)

            if chunk_data is not None:
                number_of_photons, idxs, energies = chunk_data
                photons["NumberOfPhotons"].append(number_of_photons)
                photons["Energy"].append(ds.arr(energies, "keV"))
                photons["x"].append((chunk[p_fields[0]][idxs]-parameters["center"][0]).in_units("kpc"))
                photons["y"].append((chunk[p_fields[1]][idxs]-parameters["center"][1]).in_units("kpc"))
                photons["z"].append((chunk[p_fields[2]][idxs]-parameters["center"][2]).in_units("kpc"))
                photons["vx"].append(chunk[v_fields[0]][idxs].in_units("km/s"))
                photons["vy"].append(chunk[v_fields[1]][idxs].in_units("km/s"))
                photons["vz"].append(chunk[v_fields[2]][idxs].in_units("km/s"))
                photons["dx"].append(chunk[w_field][idxs].in_units("kpc"))

        source_model.cleanup_model()

        concatenate_photons(photons)

        mylog.info("Finished generating photons.")
        mylog.info("Number of photons generated: %d" % int(np.sum(photons["NumberOfPhotons"])))
        mylog.info("Number of cells with photons: %d" % len(photons["x"]))

        return cls(photons, parameters, cosmo)

    def write_h5_file(self, photonfile):
        """
        Write the photons to the HDF5 file *photonfile*.
        """

        if parallel_capable:

            mpi_long = get_mpi_type("int64")
            mpi_double = get_mpi_type("float64")

            local_num_cells = len(self.photons["x"])
            sizes_c = comm.comm.gather(local_num_cells, root=0)

            local_num_photons = np.sum(self.photons["NumberOfPhotons"])
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
                n_ph = np.zeros(num_cells, dtype="uint64")
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
            comm.comm.Gatherv([self.photons["NumberOfPhotons"], local_num_cells, mpi_long],
                              [n_ph, (sizes_c, disps_c), mpi_long], root=0)
            comm.comm.Gatherv([self.photons["Energy"].d, local_num_photons, mpi_double],
                              [e, (sizes_p, disps_p), mpi_double], root=0)

        else:

            x = self.photons["x"].d
            y = self.photons["y"].d
            z = self.photons["z"].d
            vx = self.photons["vx"].d
            vy = self.photons["vy"].d
            vz = self.photons["vz"].d
            dx = self.photons["dx"].d
            n_ph = self.photons["NumberOfPhotons"]
            e = self.photons["Energy"].d

        if comm.rank == 0:

            f = h5py.File(photonfile, "w")

            # Parameters

            p = f.create_group("parameters")
            p.create_dataset("fid_area", data=float(self.parameters["FiducialArea"]))
            p.create_dataset("fid_exp_time", data=float(self.parameters["FiducialExposureTime"]))
            p.create_dataset("fid_redshift", data=self.parameters["FiducialRedshift"])
            p.create_dataset("hubble", data=self.parameters["HubbleConstant"])
            p.create_dataset("omega_matter", data=self.parameters["OmegaMatter"])
            p.create_dataset("omega_lambda", data=self.parameters["OmegaLambda"])
            p.create_dataset("fid_d_a", data=float(self.parameters["FiducialAngularDiameterDistance"]))
            p.create_dataset("dimension", data=self.parameters["Dimension"])
            p.create_dataset("width", data=float(self.parameters["Width"]))
            p.create_dataset("data_type", data=self.parameters["DataType"])

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

    def project_photons(self, normal, area_new=None, exp_time_new=None,
                        redshift_new=None, dist_new=None,
                        absorb_model=None, sky_center=None,
                        responses=None, convolve_energies=False, 
                        no_shifting=False, north_vector=None, 
                        prng=np.random):
        r"""
        Projects photons onto an image plane given a line of sight.

        Parameters
        ----------
        normal : character or array_like
            Normal vector to the plane of projection. If "x", "y", or "z", will
            assume to be along that axis (and will probably be faster). Otherwise,
            should be an off-axis normal vector, e.g [1.0,2.0,-3.0]
        area_new : float, (value, unit) tuple, YTQuantity, or string, optional
            New value for the effective area of the detector. A numeric value, if
            units are not specified, is assumed to be in cm**2. A string value
            indicates the name of an ARF file. If *responses* are specified the
            value of this keyword is ignored.
        exp_time_new : float, (value, unit) tuple, or YTQuantity, optional
            The new value for the exposure time. If units are not specified
            it is assumed to be in seconds.
        redshift_new : float, optional
            The new value for the cosmological redshift.
        dist_new : float, (value, unit) tuple, or YTQuantity, optional
            The new value for the angular diameter distance, used for nearby sources.
            This may be optionally supplied instead of it being determined from the
            cosmology. If units are not specified, it is assumed to be in Mpc.
        absorb_model : 'yt.analysis_modules.photon_simulator.PhotonModel`, optional
            A model for galactic absorption.
        sky_center : array_like, optional
            Center RA, Dec of the events in degrees.
        responses : list of strings, optional
            The names of the ARF and/or RMF files to convolve the photons with.
        convolve_energies : boolean, optional
            If this is set, the photon energies will be convolved with the RMF.
        no_shifting : boolean, optional
            If set, the photon energies will not be Doppler shifted.
        north_vector : a sequence of floats
            A vector defining the 'up' direction. This option sets the orientation of
            the plane of projection. If not set, an arbitrary grid-aligned north_vector
            is chosen. Ignored in the case where a particular axis (e.g., "x", "y", or
            "z") is explicitly specified.
        prng : NumPy `RandomState` object or numpy.random
            A pseudo-random number generator. Typically will only be specified if you
            have a reason to generate the same set of random numbers, such as for a
            test. Default is the numpy.random module.

        Examples
        --------
        >>> L = np.array([0.1,-0.2,0.3])
        >>> events = my_photons.project_photons(L, area_new="sim_arf.fits",
        ...                                     redshift_new=0.05,
        ...                                     psf_sigma=0.01)
        """

        if redshift_new is not None and dist_new is not None:
            mylog.error("You may specify a new redshift or distance, "+
                        "but not both!")

        if sky_center is None:
            sky_center = YTArray([30.,45.], "degree")
        else:
            sky_center = YTArray(sky_center, "degree")

        dx = self.photons["dx"].d
        nx = self.parameters["Dimension"]

        if not isinstance(normal, string_types):
            L = np.array(normal)
            orient = Orientation(L, north_vector=north_vector)
            x_hat = orient.unit_vectors[0]
            y_hat = orient.unit_vectors[1]
            z_hat = orient.unit_vectors[2]

        n_ph = self.photons["NumberOfPhotons"]
        n_ph_tot = n_ph.sum()

        eff_area = None

        parameters = {}

        if responses is not None:
            responses = ensure_list(responses)
            parameters["ARF"] = responses[0]
            if len(responses) == 2:
                parameters["RMF"] = responses[1]
            area_new = parameters["ARF"]

        zobs0 = self.parameters["FiducialRedshift"]
        D_A0 = self.parameters["FiducialAngularDiameterDistance"]
        scale_factor = 1.0

        # If we use an RMF, figure out where the response matrix actually is.
        if "RMF" in parameters:
            rmf = _astropy.pyfits.open(parameters["RMF"])
            if "MATRIX" in rmf:
                mat_key = "MATRIX"
            elif "SPECRESP MATRIX" in rmf:
                mat_key = "SPECRESP MATRIX"
            else:
                raise RuntimeError("Cannot find the response matrix in the RMF "
                                   "file %s! " % parameters["RMF"]+"It should "
                                   "be named \"MATRIX\" or \"SPECRESP MATRIX\".")
            rmf.close()
        else:
            mat_key = None

        if (exp_time_new is None and area_new is None and
            redshift_new is None and dist_new is None):
            my_n_obs = n_ph_tot
            zobs = zobs0
            D_A = D_A0
        else:
            if exp_time_new is None:
                Tratio = 1.
            else:
                Tratio = parse_value(exp_time_new, "s")/self.parameters["FiducialExposureTime"]
            if area_new is None:
                Aratio = 1.
            elif isinstance(area_new, string_types):
                if comm.rank == 0:
                    mylog.info("Using energy-dependent effective area: %s" % (parameters["ARF"]))
                f = _astropy.pyfits.open(area_new)
                earf = 0.5*(f["SPECRESP"].data.field("ENERG_LO")+f["SPECRESP"].data.field("ENERG_HI"))
                eff_area = np.nan_to_num(f["SPECRESP"].data.field("SPECRESP"))
                if "RMF" in parameters:
                    weights = self._normalize_arf(parameters["RMF"], mat_key)
                    eff_area *= weights
                else:
                    mylog.warning("You specified an ARF but not an RMF. This is ok if the "+
                                  "responses are normalized properly. If not, you may "+
                                  "get inconsistent results.")
                f.close()
                Aratio = eff_area.max()/self.parameters["FiducialArea"].v
            else:
                mylog.info("Using constant effective area.")
                Aratio = parse_value(area_new, "cm**2")/self.parameters["FiducialArea"]
            if redshift_new is None and dist_new is None:
                Dratio = 1.
                zobs = zobs0
                D_A = D_A0
            else:
                if redshift_new is None:
                    zobs = 0.0
                    D_A = parse_value(dist_new, "Mpc")
                else:
                    zobs = redshift_new
                    D_A = self.cosmo.angular_diameter_distance(0.0,zobs).in_units("Mpc")
                    scale_factor = (1.+zobs0)/(1.+zobs)
                Dratio = D_A0*D_A0*(1.+zobs0)**3 / \
                         (D_A*D_A*(1.+zobs)**3)
            fak = Aratio*Tratio*Dratio
            if fak > 1:
                raise ValueError("This combination of requested parameters results in "
                                 "%g%% more photons collected than are " % (100.*(fak-1.)) +
                                 "available in the sample. Please reduce the collecting "
                                 "area, exposure time, or increase the distance/redshift "
                                 "of the object. Alternatively, generate a larger sample "
                                 "of photons.")
            my_n_obs = np.uint64(n_ph_tot*fak)

        n_obs_all = comm.mpi_allreduce(my_n_obs)
        if comm.rank == 0:
            mylog.info("Total number of photons to use: %d" % n_obs_all)

        if my_n_obs == n_ph_tot:
            idxs = np.arange(my_n_obs, dtype='uint64')
        else:
            idxs = prng.permutation(n_ph_tot)[:my_n_obs].astype("uint64")
        obs_cells = np.searchsorted(self.p_bins, idxs, side='right')-1
        delta = dx[obs_cells]

        if isinstance(normal, string_types):

            if self.parameters["DataType"] == "cells":
                xsky = prng.uniform(low=-0.5, high=0.5, size=my_n_obs)
                ysky = prng.uniform(low=-0.5, high=0.5, size=my_n_obs)
            elif self.parameters["DataType"] == "particles":
                xsky = prng.normal(loc=0.0, scale=1.0, size=my_n_obs)
                ysky = prng.normal(loc=0.0, scale=1.0, size=my_n_obs)
            xsky *= delta
            ysky *= delta
            xsky += self.photons[axes_lookup[normal][0]].d[obs_cells]
            ysky += self.photons[axes_lookup[normal][1]].d[obs_cells]

            if not no_shifting:
                vz = self.photons["v%s" % normal]

        else:

            if self.parameters["DataType"] == "cells":
                x = prng.uniform(low=-0.5, high=0.5, size=my_n_obs)
                y = prng.uniform(low=-0.5, high=0.5, size=my_n_obs)
                z = prng.uniform(low=-0.5, high=0.5, size=my_n_obs)
            elif self.parameters["DataType"] == "particles":
                x = prng.normal(loc=0.0, scale=1.0, size=my_n_obs)
                y = prng.normal(loc=0.0, scale=1.0, size=my_n_obs)
                z = prng.normal(loc=0.0, scale=1.0, size=my_n_obs)

            if not no_shifting:
                vz = self.photons["vx"]*z_hat[0] + \
                     self.photons["vy"]*z_hat[1] + \
                     self.photons["vz"]*z_hat[2]

            x *= delta
            y *= delta
            z *= delta
            x += self.photons["x"].d[obs_cells]
            y += self.photons["y"].d[obs_cells]
            z += self.photons["z"].d[obs_cells]

            xsky = x*x_hat[0] + y*x_hat[1] + z*x_hat[2]
            ysky = x*y_hat[0] + y*y_hat[1] + z*y_hat[2]

        if no_shifting:
            eobs = self.photons["Energy"][idxs]
        else:
            shift = -vz.in_cgs()/clight
            shift = np.sqrt((1.-shift)/(1.+shift))
            eobs = self.photons["Energy"][idxs]*shift[obs_cells]
        eobs *= scale_factor

        if absorb_model is None:
            not_abs = np.ones(eobs.shape, dtype='bool')
        else:
            mylog.info("Absorbing.")
            absorb_model.prepare_spectrum()
            emid = absorb_model.emid
            aspec = absorb_model.get_spectrum()
            absorb = np.interp(eobs, emid, aspec, left=0.0, right=0.0)
            randvec = aspec.max()*prng.uniform(size=eobs.shape)
            not_abs = randvec < absorb
            absorb_model.cleanup_spectrum()

        if eff_area is None:
            detected = np.ones(eobs.shape, dtype='bool')
        else:
            mylog.info("Applying energy-dependent effective area.")
            earea = np.interp(eobs, earf, eff_area, left=0.0, right=0.0)
            randvec = eff_area.max()*prng.uniform(size=eobs.shape)
            detected = randvec < earea

        detected = np.logical_and(not_abs, detected)

        events = {}

        dx_min = self.parameters["Width"]/self.parameters["Dimension"]
        dtheta = YTQuantity(np.rad2deg(dx_min/D_A), "degree")

        events["xpix"] = xsky[detected]/dx_min.v + 0.5*(nx+1)
        events["ypix"] = ysky[detected]/dx_min.v + 0.5*(nx+1)
        events["eobs"] = eobs[detected]

        events = comm.par_combine_object(events, datatype="dict", op="cat")

        num_events = len(events["xpix"])

        if comm.rank == 0:
            mylog.info("Total number of observed photons: %d" % num_events)

        if "RMF" in parameters and convolve_energies:
            events, info = self._convolve_with_rmf(parameters["RMF"], events, 
                                                   mat_key, prng)
            parameters.update(info)

        if exp_time_new is None:
            parameters["ExposureTime"] = self.parameters["FiducialExposureTime"]
        else:
            parameters["ExposureTime"] = exp_time_new
        if area_new is None:
            parameters["Area"] = self.parameters["FiducialArea"]
        else:
            parameters["Area"] = area_new
        parameters["Redshift"] = zobs
        parameters["AngularDiameterDistance"] = D_A.in_units("Mpc")
        parameters["sky_center"] = sky_center
        parameters["pix_center"] = np.array([0.5*(nx+1)]*2)
        parameters["dtheta"] = dtheta

        return EventList(events, parameters)

    def _normalize_arf(self, respfile, mat_key):
        rmf = _astropy.pyfits.open(respfile)
        table = rmf[mat_key]
        weights = np.array([w.sum() for w in table.data["MATRIX"]])
        rmf.close()
        return weights

    def _convolve_with_rmf(self, respfile, events, mat_key, prng):
        """
        Convolve the events with a RMF file.
        """
        mylog.info("Reading response matrix file (RMF): %s" % (respfile))

        hdulist = _astropy.pyfits.open(respfile)

        tblhdu = hdulist[mat_key]
        n_de = len(tblhdu.data["ENERG_LO"])
        mylog.info("Number of energy bins in RMF: %d" % (n_de))
        mylog.info("Energy limits: %g %g" % (min(tblhdu.data["ENERG_LO"]),
                                             max(tblhdu.data["ENERG_HI"])))

        tblhdu2 = hdulist["EBOUNDS"]
        n_ch = len(tblhdu2.data["CHANNEL"])
        mylog.info("Number of channels in RMF: %d" % (n_ch))

        eidxs = np.argsort(events["eobs"])

        phEE = events["eobs"][eidxs].d

        detectedChannels = []

        # run through all photon energies and find which bin they go in
        k = 0
        fcurr = 0
        last = len(phEE)

        pbar = get_pbar("Scattering energies with RMF", last)

        for low,high in zip(tblhdu.data["ENERG_LO"],tblhdu.data["ENERG_HI"]):
            # weight function for probabilities from RMF
            weights = np.nan_to_num(np.float64(tblhdu.data[k]["MATRIX"][:]))
            weights /= weights.sum()
            # build channel number list associated to array value,
            # there are groups of channels in rmfs with nonzero probabilities
            trueChannel = []
            f_chan = np.nan_to_num(tblhdu.data["F_CHAN"][k])
            n_chan = np.nan_to_num(tblhdu.data["N_CHAN"][k])
            n_grp = np.nan_to_num(tblhdu.data["N_CHAN"][k])
            if not iterable(f_chan):
                f_chan = [f_chan]
                n_chan = [n_chan]
                n_grp  = [n_grp]
            for start,nchan in zip(f_chan, n_chan):
                end = start + nchan
                if start == end:
                    trueChannel.append(start)
                else:
                    for j in range(start,end):
                        trueChannel.append(j)
            if len(trueChannel) > 0:
                for q in range(fcurr,last):
                    if phEE[q] >= low and phEE[q] < high:
                        channelInd = prng.choice(len(weights), p=weights)
                        fcurr += 1
                        pbar.update(fcurr)
                        detectedChannels.append(trueChannel[channelInd])
                    if phEE[q] >= high:
                        break
            k += 1
        pbar.finish()

        dchannel = np.array(detectedChannels)

        events["xpix"] = events["xpix"][eidxs]
        events["ypix"] = events["ypix"][eidxs]
        events["eobs"] = YTArray(phEE, "keV")
        events[tblhdu.header["CHANTYPE"]] = dchannel.astype(int)

        info = {"ChannelType" : tblhdu.header["CHANTYPE"],
                "Telescope" : tblhdu.header["TELESCOP"],
                "Instrument" : tblhdu.header["INSTRUME"]}

        info["Mission"] = tblhdu.header.get("MISSION","")

        return events, info