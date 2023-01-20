"""
Classes for generating lists of photons
"""
import h5py
import numpy as np
from soxs import __version__ as soxs_version
from soxs.utils import parse_prng
from tqdm.auto import tqdm
from unyt.array import unyt_array
from yt import __version__ as yt_version
from yt.utilities.cosmology import Cosmology
from yt.utilities.orientation import Orientation
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    communication_system,
    parallel_objects,
)
from yt.utilities.physical_constants import clight

from pyxsim import __version__ as pyxsim_version
from pyxsim.lib.sky_functions import (
    doppler_shift,
    pixel_to_cel,
    scatter_events,
    scatter_events_allsky,
)
from pyxsim.spectral_models import absorb_models
from pyxsim.utils import mylog, parse_value

comm = communication_system.communicators[-1]

init_chunk = 100000


def determine_fields(ds, source_type, point_sources):
    from yt.geometry.particle_geometry_handler import ParticleIndex

    # Figure out if this is a particle field or otherwise
    ptype = (
        (source_type in ds.particle_types)
        | (source_type in ds.known_filters)
        | ((source_type == "gas") & isinstance(ds.index, ParticleIndex))
    )
    if ptype:
        ppos = [f"particle_position_{ax}" for ax in "xyz"]
        pvel = [f"particle_velocity_{ax}" for ax in "xyz"]
        if source_type in ds.known_filters:
            if ds.known_filters[source_type].filtered_type == "gas":
                ppos = ["x", "y", "z"]
                pvel = [f"velocity_{ax}" for ax in "xyz"]
        elif source_type == "gas":
            source_type = ds._sph_ptypes[0]
        position_fields = [(source_type, ppos[i]) for i in range(3)]
        velocity_fields = [(source_type, pvel[i]) for i in range(3)]
        if source_type in ds._sph_ptypes:
            width_field = (source_type, "smoothing_length")
        else:
            width_field = None
    else:
        position_fields = [("index", ax) for ax in "xyz"]
        velocity_fields = [(source_type, f"velocity_{ax}") for ax in "xyz"]
        width_field = ("index", "dx")
    if point_sources:
        width_field = None
    return position_fields, velocity_fields, width_field


def find_object_bounds(data_source):
    """
    This logic is required to determine the bounds of the object, which is
    solely for fixing coordinates at periodic boundaries
    """
    if hasattr(data_source, "base_object"):
        # This is a cut region so we'll figure out
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
        mylog.warning(
            "You are using a region that is not currently "
            "supported for straddling periodic boundaries. "
            "Check to make sure that your region does not "
            "do so."
        )
        le = data_source.ds.domain_left_edge
        re = data_source.ds.domain_right_edge

    return le.to_value("kpc"), re.to_value("kpc")


def make_photons(
    photon_prefix,
    data_source,
    redshift,
    area,
    exp_time,
    source_model,
    point_sources=False,
    parameters=None,
    center=None,
    dist=None,
    cosmology=None,
    velocity_fields=None,
    bulk_velocity=None,
    observer="external",
    fields_to_keep=None,
):
    r"""
    Write a photon list dataset to disk from a yt data source and assuming a
    source model for the photons. The redshift, collecting area, exposure time,
    and cosmology are stored in the *parameters* dictionary which is passed to
    the *source_model* function.

    Parameters
    ----------
    photon_prefix : string
        The prefix of the filename(s) to contain the photon list. If run in
        serial, the filename will be "{photon_prefix}.h5", if run in
        parallel, the filenames will be "{photon_prefix}.{mpi_rank}.h5".
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
    source_model : :class:`~pyxsim.source_models.sources.SourceModel`
        A source model used to generate the photons.
    point_sources : boolean, optional
        If True, the photons will be assumed to be generated from the exact
        positions of the cells or particles and not smeared around within
        a volume. Default: False
    parameters : dict, optional
        A dictionary of parameters to be passed for the source model to use,
        if necessary.
    center : string or array_like, optional
        The origin of the photon spatial coordinates. Accepts "c", "max", or
        a coordinate. If array-like and without units, it is assumed to be in
        units of kpc. If not specified, pyxsim attempts to use the "center"
        field parameter of the data_source.
    dist : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`, optional
        The angular diameter distance, used for nearby sources. This may be
        optionally supplied instead of it being determined from the
        *redshift* and given *cosmology*. If units are not specified, it is
        assumed to be in kpc. To use this, the redshift must be set to zero.
    cosmology : :class:`~yt.utilities.cosmology.Cosmology`, optional
        Cosmological information. If not supplied, we try to get the
        cosmology from the dataset. Otherwise, LCDM with the default yt
        parameters is assumed.
    velocity_fields : list of fields
        The yt fields to use for the velocity. If not specified, the
        following will be assumed:
        ['velocity_x', 'velocity_y', 'velocity_z'] for grid datasets
        ['particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z'] for particle datasets
    bulk_velocity : array-like, optional
        A 3-element array or list specifying the local velocity frame of
        reference. If not a :class:`~yt.units.yt_array.YTArray`, it is assumed
        to have units of km/s. Default: [0.0, 0.0, 0.0] km/s.

    Returns
    -------
    A tuple of two integers, the number of photons and the number of cells with
    photons

    Examples
    --------
    >>> thermal_model = pyxsim.CIESourceModel(0.1, 10.0, 1000, 0.3)
    >>> redshift = 0.05
    >>> area = 6000.0 # assumed here in cm**2
    >>> time = 2.0e5 # assumed here in seconds
    >>> sp = ds.sphere("c", (500., "kpc"))
    >>> n_photons, n_cells = pyxsim.make_photons(sp, redshift, area,
    ...                                          time, thermal_model)
    """
    ds = data_source.ds

    if photon_prefix.endswith(".h5"):
        photon_prefix = photon_prefix[:-3]

    if comm.size > 1:
        photon_file = f"{photon_prefix}.{comm.rank:04d}.h5"
    else:
        photon_file = f"{photon_prefix}.h5"

    if parameters is None:
        parameters = {}
    if cosmology is None:
        if hasattr(ds, "cosmology"):
            cosmo = ds.cosmology
        else:
            cosmo = Cosmology()
    else:
        cosmo = cosmology

    if observer == "external":
        if dist is None:
            if redshift <= 0.0:
                msg = (
                    "If redshift <= 0.0, you must specify a distance to the "
                    "source using the 'dist' argument!"
                )
                mylog.error(msg)
                raise ValueError(msg)
            D_A = cosmo.angular_diameter_distance(0.0, redshift).to("Mpc")
        else:
            D_A = parse_value(dist, "kpc")
            if redshift > 0.0:
                mylog.warning(
                    "Redshift must be zero for nearby sources. "
                    "Resetting redshift to 0.0."
                )
                redshift = 0.0
    else:
        D_A = parse_value(0.0, "kpc")
        if redshift > 0.0:
            mylog.warning(
                "Redshift must be zero for internal observers. "
                "Resetting redshift to 0.0."
            )
            redshift = 0.0

    if isinstance(center, str):
        if center == "center" or center == "c":
            parameters["center"] = ds.domain_center
        elif center == "max" or center == "m":
            parameters["center"] = ds.find_max("density")[-1]
    elif isinstance(center, unyt_array):
        parameters["center"] = center.in_units("code_length")
    elif isinstance(center, tuple):
        if center[0] == "min":
            parameters["center"] = ds.find_min(center[1])[-1]
        elif center[0] == "max":
            parameters["center"] = ds.find_max(center[1])[-1]
        else:
            raise RuntimeError
    elif isinstance(center, (list, np.ndarray)):
        parameters["center"] = ds.arr(center, "code_length")
    elif center is None:
        if hasattr(data_source, "left_edge"):
            parameters["center"] = 0.5 * (
                data_source.left_edge + data_source.right_edge
            )
        else:
            parameters["center"] = data_source.get_field_parameter("center")

    if bulk_velocity is None:
        bulk_velocity = ds.arr([0.0] * 3, "km/s")
    elif isinstance(bulk_velocity, unyt_array):
        bulk_velocity = bulk_velocity.to("km/s")
    elif isinstance(bulk_velocity, (list, np.ndarray)):
        bulk_velocity = ds.arr(bulk_velocity, "km/s")
    parameters["bulk_velocity"] = bulk_velocity

    parameters["fid_exp_time"] = parse_value(exp_time, "s")
    parameters["fid_area"] = parse_value(area, "cm**2")
    parameters["fid_redshift"] = redshift
    parameters["observer"] = observer
    parameters["hubble"] = cosmo.hubble_constant
    parameters["omega_matter"] = cosmo.omega_matter
    parameters["omega_lambda"] = cosmo.omega_lambda
    parameters["center"].convert_to_units("kpc")
    parameters["fid_d_a"] = D_A

    if observer == "external":
        if redshift > 0.0:
            mylog.info(
                "Cosmology: h = %g, omega_matter = %g, omega_lambda = %g",
                cosmo.hubble_constant,
                cosmo.omega_matter,
                cosmo.omega_lambda,
            )
        else:
            mylog.info("Observing local source at distance %g.", D_A)
    else:
        mylog.info("The observer is internal to the source.")

    local_exp_time = parameters["fid_exp_time"].v
    D_A = parameters["fid_d_a"].to_value("cm")
    dist_fac = 1.0 / (4.0 * np.pi)
    if observer == "external":
        dist_fac /= D_A * D_A * (1.0 + redshift) ** 2

    spectral_norm = parameters["fid_area"].v * local_exp_time * dist_fac

    dw = ds.domain_width.to_value("kpc")
    le, re = find_object_bounds(data_source)
    c = parameters["center"].to_value("kpc")

    source_model.setup_model("photons", data_source, redshift)

    p_fields, v_fields, w_field = determine_fields(
        ds, source_model.ftype, point_sources
    )

    if velocity_fields is not None:
        v_fields = velocity_fields

    fields_store = []
    if fields_to_keep is not None:
        for field in fields_to_keep:
            fd = ds._get_field_info(field)
            fields_store.append(fd.name)

    if p_fields[0] == ("index", "x"):
        parameters["data_type"] = "cells"
    else:
        parameters["data_type"] = "particles"

    source_model.set_pv(p_fields, v_fields, le, re, dw, c, ds.periodicity, observer)

    f = h5py.File(photon_file, "w")

    # Info
    info = f.create_group("info")
    info.attrs["yt_version"] = yt_version
    info.attrs["pyxsim_version"] = pyxsim_version
    info.attrs["soxs_version"] = soxs_version
    info.attrs["dataset"] = str(data_source.ds)
    info.attrs["data_source"] = str(data_source)
    info.attrs["source_model"] = repr(source_model)

    # Parameters

    p = f.create_group("parameters")
    p.create_dataset("fid_area", data=float(parameters["fid_area"]))
    p.create_dataset("fid_exp_time", data=float(parameters["fid_exp_time"]))
    p.create_dataset("fid_redshift", data=parameters["fid_redshift"])
    p.create_dataset("hubble", data=parameters["hubble"])
    p.create_dataset("omega_matter", data=parameters["omega_matter"])
    p.create_dataset("omega_lambda", data=parameters["omega_lambda"])
    p.create_dataset("fid_d_a", data=float(parameters["fid_d_a"]))
    p.create_dataset("data_type", data=parameters["data_type"])
    p.create_dataset("observer", data=parameters["observer"])
    p.create_dataset("center", data=parameters["center"].d)
    p.create_dataset("bulk_velocity", data=parameters["bulk_velocity"].d)
    p.create_dataset("velocity_fields", data=np.array(v_fields).astype("S"))

    n_cells = 0
    n_photons = 0
    c_offset = 0
    p_offset = 0
    c_size = init_chunk
    p_size = init_chunk

    cell_fields = ["x", "y", "z", "vx", "vy", "vz", "num_photons", "dx"]
    if len(fields_store) > 0:
        for field in fields_store:
            cell_fields.append(field[1])

    d = f.create_group("data")
    for field in cell_fields + ["energy"]:
        if field == "num_photons":
            dtype = "int64"
        else:
            dtype = "float64"
        d.create_dataset(
            field,
            data=np.zeros(init_chunk, dtype=dtype),
            maxshape=(None,),
            dtype=dtype,
            chunks=True,
        )

    f.flush()

    for chunk in parallel_objects(data_source.chunks([], "io")):

        chunk_data = source_model.process_data("photons", chunk, spectral_norm)

        if chunk_data is not None:

            chunk_nc, number_of_photons, idxs, energies = chunk_data

            chunk_nph = np.sum(number_of_photons)

            if chunk_nph == 0:
                continue

            if c_size < n_cells + chunk_nc:
                while chunk_nc + n_cells > c_size:
                    c_size *= 2
                for field in cell_fields:
                    d[field].resize((c_size,))

            if p_size < n_photons + chunk_nph:
                while chunk_nph + n_photons > p_size:
                    p_size *= 2
                d["energy"].resize((p_size,))

            for i, ax in enumerate("xyz"):
                pos = chunk[p_fields[i]][idxs].to_value("kpc")
                # Fix photon coordinates for regions crossing a periodic boundary
                if ds.periodicity[i]:
                    tfl = pos < le[i]
                    tfr = pos > re[i]
                    pos[tfl] += dw[i]
                    pos[tfr] -= dw[i]

                vel = chunk[v_fields[i]][idxs].to_value("km/s")
                # Coordinates are centered
                d[ax][c_offset : c_offset + chunk_nc] = pos - c[i]
                # Velocities have the bulk velocity subtracted off
                d[f"v{ax}"][c_offset : c_offset + chunk_nc] = vel - bulk_velocity.v[i]

            d["num_photons"][c_offset : c_offset + chunk_nc] = number_of_photons
            d["energy"][p_offset : p_offset + chunk_nph] = energies

            if w_field is None:
                d["dx"][c_offset : c_offset + chunk_nc] = 0.0
            else:
                d["dx"][c_offset : c_offset + chunk_nc] = chunk[w_field][idxs].to_value(
                    "kpc"
                )

            for field in fields_store:
                d[field[1]][c_offset : c_offset + chunk_nc] = chunk[field][idxs]

            n_cells += chunk_nc
            n_photons += chunk_nph
            c_offset = n_cells
            p_offset = n_photons

        f.flush()

    if c_size > n_cells:
        for field in cell_fields:
            d[field].resize((n_cells,))

    if p_size > n_photons:
        d["energy"].resize((n_photons,))

    f.close()

    source_model.cleanup_model("photons")

    all_nphotons = comm.mpi_allreduce(n_photons)
    all_ncells = comm.mpi_allreduce(n_cells)

    mylog.info("Finished generating photons.")
    mylog.info("Number of photons generated: %d", all_nphotons)
    mylog.info("Number of cells with photons: %d", all_ncells)

    return all_nphotons, all_ncells


def _project_photons(
    obs,
    photon_prefix,
    event_prefix,
    normal,
    sky_center,
    absorb_model=None,
    nH=None,
    abund_table="angr",
    no_shifting=False,
    north_vector=None,
    flat_sky=False,
    sigma_pos=None,
    kernel="top_hat",
    save_los=False,
    prng=None,
):

    from yt.funcs import ensure_numpy_array

    prng = parse_prng(prng)

    if photon_prefix.endswith(".h5"):
        photon_prefix = photon_prefix[:-3]

    if event_prefix.endswith(".h5"):
        event_prefix = event_prefix[:-3]

    if not isinstance(normal, str):
        L = np.array(normal)
        orient = Orientation(L, north_vector=north_vector)
        x_hat = orient.unit_vectors[0]
        y_hat = orient.unit_vectors[1]
        z_hat = orient.unit_vectors[2]
        north_vector = orient.north_vector
    else:
        x_hat = np.zeros(3)
        y_hat = np.zeros(3)
        z_hat = np.zeros(3)
        north_vector = None

    if comm.size > 1:
        photon_file = f"{photon_prefix}.{comm.rank:04d}.h5"
        event_file = f"{event_prefix}.{comm.rank:04d}.h5"
    else:
        photon_file = f"{photon_prefix}.h5"
        event_file = f"{event_prefix}.h5"

    sky_center = ensure_numpy_array(sky_center)

    scale_shift = -1.0 / clight.to_value("km/s")
    scale_shift2 = scale_shift * scale_shift

    if isinstance(absorb_model, str):
        if absorb_model not in absorb_models:
            raise KeyError(f"{absorb_model} is not a known absorption model!")
        absorb_model = absorb_models[absorb_model]
    if absorb_model is not None:
        if nH is None:
            raise RuntimeError(
                "You specified an absorption model, but didn't "
                "specify a value for nH!"
            )
        absorb_model = absorb_model(nH, abund_table=abund_table)
        if comm.rank == 0:
            mylog.info(
                "Foreground galactic absorption: using the %s model and nH = %g.",
                absorb_model._name,
                nH,
            )
    abs_model_name = absorb_model._name if absorb_model else "none"
    if nH is None:
        nH = 0.0

    f = h5py.File(photon_file, "r")

    p = f["parameters"]

    data_type = p["data_type"].asstr()[()]
    if "observer" in p:
        observer = p["observer"].asstr()[()]
    else:
        observer = "external"
    if obs != observer:
        which_func = {"external": "", "internal": "_allsky"}
        raise RuntimeError(
            f"The function 'project_photons{which_func['obs']}' "
            f"does not work with '{observer}' photon lists!"
        )

    if observer == "internal" and isinstance(normal, str):
        raise RuntimeError(
            "Must specify a vector for 'normal' if you are "
            "doing an 'internal' observation!"
        )

    if sigma_pos is not None and data_type == "particles":
        raise RuntimeError(
            "The 'smooth_positions' argument should "
            "not be used with particle-based datasets!"
        )

    d = f["data"]

    D_A = p["fid_d_a"][()] * 1.0e3

    if d["energy"].size == 0:

        mylog.warning("No photons are in file %s, so I am done.", photon_file)
        n_events = 0

    else:

        fe = h5py.File(event_file, "w")

        ie = fe.create_group("info")
        ie.attrs["pyxsim_version"] = pyxsim_version
        ie.attrs["yt_version"] = yt_version
        ie.attrs["soxs_version"] = soxs_version
        ie.attrs["photon_file"] = photon_file

        pe = fe.create_group("parameters")
        pe.create_dataset("exp_time", data=float(p["fid_exp_time"][()]))
        pe.create_dataset("area", data=float(p["fid_area"][()]))
        pe.create_dataset("sky_center", data=sky_center)
        pe.create_dataset("observer", data=observer)
        pe.create_dataset("no_shifting", data=int(no_shifting))
        pe.create_dataset("flat_sky", data=int(flat_sky))
        pe.create_dataset("normal", data=normal)
        if north_vector is not None:
            pe.create_dataset("north_vector", data=north_vector)
        pe.create_dataset("absoption_model", data=abs_model_name)
        if absorb_model is not None:
            pe.create_dataset("nH", data=nH)
            pe.create_dataset("abund_table", data=abund_table)
        if sigma_pos is not None:
            pe.create_dataset("sigma_pos", data=sigma_pos)
        pe.create_dataset("kernel", data=kernel)
        event_fields = ["xsky", "ysky", "eobs"]
        if save_los:
            event_fields.append("los")

        n_events = 0
        e_offset = 0
        e_size = init_chunk
        cell_chunk = init_chunk
        start_e = 0

        de = fe.create_group("data")
        for field in event_fields:
            de.create_dataset(
                field, data=np.zeros(init_chunk), maxshape=(None,), chunks=True
            )

        if isinstance(normal, str):
            norm = "xyz".index(normal)
        else:
            norm = normal

        n_cells = d["num_photons"].size

        pbar = tqdm(
            leave=True, total=n_cells, desc="Projecting photons from cells/particles "
        )

        for start_c in range(0, n_cells, cell_chunk):

            end_c = min(start_c + cell_chunk, n_cells)

            n_ph = d["num_photons"][start_c:end_c]
            x = d["x"][start_c:end_c]
            y = d["y"][start_c:end_c]
            z = d["z"][start_c:end_c]
            dx = d["dx"][start_c:end_c]
            end_e = start_e + n_ph.sum()
            eobs = d["energy"][start_e:end_e]

            if observer == "internal":
                r = np.sqrt(x * x + y * y + z * z)
            else:
                r = None

            if not no_shifting:
                if observer == "internal":
                    vn = (
                        -(
                            d["vx"][start_c:end_c] * x
                            + d["vy"][start_c:end_c] * y
                            + d["vz"][start_c:end_c] * z
                        )
                        / r
                    )
                else:
                    if isinstance(normal, str):
                        vn = d[f"v{normal}"][start_c:end_c]
                    else:
                        vn = (
                            d["vx"][start_c:end_c] * z_hat[0]
                            + d["vy"][start_c:end_c] * z_hat[1]
                            + d["vz"][start_c:end_c] * z_hat[2]
                        )
                v2 = (
                    d["vx"][start_c:end_c] * d["vx"][start_c:end_c]
                    + d["vy"][start_c:end_c] * d["vy"][start_c:end_c]
                    + d["vz"][start_c:end_c] * d["vz"][start_c:end_c]
                )
                doppler_shift(vn * scale_shift, v2 * scale_shift2, n_ph, eobs)

            if absorb_model is None:
                det = np.ones(eobs.size, dtype="bool")
                num_det = eobs.size
            else:
                det = absorb_model.absorb_photons(eobs, prng=prng)
                num_det = det.sum()

            if num_det > 0:

                if observer == "external":

                    if data_type == "particles":
                        dx *= 0.5

                    xsky, ysky, los = scatter_events(
                        norm,
                        prng,
                        kernel,
                        data_type,
                        num_det,
                        det,
                        n_ph,
                        x,
                        y,
                        z,
                        dx,
                        x_hat,
                        y_hat,
                        z_hat,
                    )

                    if data_type == "cells" and sigma_pos is not None:
                        sigma = sigma_pos * np.repeat(dx, n_ph)[det]
                        xsky += sigma * prng.normal(loc=0.0, scale=1.0, size=num_det)
                        ysky += sigma * prng.normal(loc=0.0, scale=1.0, size=num_det)

                    xsky /= D_A
                    ysky /= D_A

                    if flat_sky:
                        xsky = sky_center[0] - np.rad2deg(xsky)
                        ysky = sky_center[1] + np.rad2deg(ysky)
                    else:
                        pixel_to_cel(xsky, ysky, sky_center)

                elif observer == "internal":

                    xsky, ysky, los = scatter_events_allsky(
                        data_type,
                        kernel,
                        prng,
                        num_det,
                        det,
                        n_ph,
                        x,
                        y,
                        z,
                        dx,
                        x_hat,
                        y_hat,
                        z_hat,
                    )

                if e_size < n_events + num_det:
                    while n_events + num_det > e_size:
                        e_size *= 2
                    for field in event_fields:
                        de[field].resize((e_size,))

                de["xsky"][e_offset : e_offset + num_det] = xsky
                de["ysky"][e_offset : e_offset + num_det] = ysky
                de["eobs"][e_offset : e_offset + num_det] = eobs[det]
                if save_los:
                    de["los"][e_offset : e_offset + num_det] = los

                n_events += num_det
                e_offset = n_events

                f.flush()

            pbar.update(end_c - start_c + 1)

            start_e = end_e

        pbar.close()

        if e_size > n_events:
            for field in event_fields:
                de[field].resize((n_events,))

        fe.close()

    f.close()

    all_nevents = comm.mpi_allreduce(n_events)

    mylog.info("Detected %d events.", all_nevents)

    return all_nevents


def project_photons(
    photon_prefix,
    event_prefix,
    normal,
    sky_center,
    absorb_model=None,
    nH=None,
    abund_table="angr",
    no_shifting=False,
    north_vector=None,
    sigma_pos=None,
    flat_sky=False,
    kernel="top_hat",
    save_los=False,
    prng=None,
):
    r"""
    Projects photons onto an image plane given a line of sight, and
    stores them in an HDF5 dataset which contains an event list.

    Parameters
    ----------
    photon_prefix : string
        The prefix of the filename(s) containing the photon list. If run in
        serial, the filename will be "{photon_prefix}.h5", if run in
        parallel, the filenames will be "{photon_prefix}.{mpi_rank}.h5".
    event_prefix : string
        The prefix of the filename(s) which will be written to contain the
        event list. If run in serial, the filename will be "{event_prefix}.h5",
        if run in parallel, the filename will be "{event_prefix}.{mpi_rank}.h5".
    normal : character or array-like
        Normal vector to the plane of projection. If "x", "y", or "z", will
        assume to be along that axis (and will probably be faster). Otherwise,
        should be an off-axis normal vector, e.g [1.0, 2.0, -3.0]
    sky_center : array-like
        Center RA, Dec of the events in degrees.
    absorb_model : string
        A model for foreground galactic absorption, to simulate the
        absorption of events before being detected. Known options for
        are "wabs" and "tbabs".
    nH : float, optional
        The foreground column density in units of 10^22 cm^{-2}. Only used
        if absorption is applied.
    abund_table : string
        The abundance table to be used for abundances in the
        absorption model (only used for TBabs). Default is set in the SOXS
        configuration file, the default for which is "angr".
        Built-in options are:
        "angr" : from Anders E. & Grevesse N. (1989, Geochimica et
        Cosmochimica Acta 53, 197)
        "aspl" : from Asplund M., Grevesse N., Sauval A.J. & Scott
        P. (2009, ARAA, 47, 481)
        "feld" : from Feldman U. (1992, Physica Scripta, 46, 202)
        "wilm" : from Wilms, Allen & McCray (2000, ApJ 542, 914
        except for elements not listed which are given zero abundance)
        "lodd" : from Lodders, K (2003, ApJ 591, 1220)
        "cl17.03" : the default abundance table in Cloudy 17.03
    no_shifting : boolean, optional
        If set, the photon energies will not be Doppler shifted. Default: False
    north_vector : a sequence of floats
        A vector defining the "up" direction. This option sets the
        orientation of the plane of projection. If not set, an arbitrary
        grid-aligned north_vector perpendicular to the normal is chosen.
        Ignored in the case where a particular axis (e.g., "x", "y", or
        "z") is explicitly specified.
    sigma_pos : float, optional
        Apply a gaussian smoothing operation to the sky positions of the
        events. This may be useful when the binned events appear blocky due
        to their uniform distribution within simulation cells. However, this
        will move the events away from their originating position on the
        sky, and so may distort surface brightness profiles and/or spectra.
        Should probably only be used for visualization purposes. Supply a
        float here to smooth with a standard deviation with this fraction
        of the cell size. Default: None
    flat_sky : boolean, optional
        If set, we assume that the sky is "flat" and RA, Dec positions are
        computed using simple linear offsets
    kernel : string, optional
        The kernel used when smoothing positions of X-rays originating from
        SPH particles, "gaussian" or "top_hat". Default: "top_hat".
    save_los : boolean, optional
        If True, save the line-of-sight positions along the projection axis in
        units of kpc to the events list. Default: False
    prng : integer or :class:`~numpy.random.RandomState` object
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is to use the :mod:`numpy.random`
        module.

    Returns
    -------
    A integer for the number of events created

    Examples
    --------
    >>> L = np.array([0.1,-0.2,0.3])
    >>> n_events = pyxsim.project_photons("my_photons.h5", "my_events.h5", L,
    ...                                   [30., 45.], absorb_model='tbabs',
    ...                                   nH=0.04)
    """
    return _project_photons(
        "external",
        photon_prefix,
        event_prefix,
        normal,
        sky_center,
        absorb_model=absorb_model,
        nH=nH,
        abund_table=abund_table,
        no_shifting=no_shifting,
        north_vector=north_vector,
        sigma_pos=sigma_pos,
        flat_sky=flat_sky,
        kernel=kernel,
        save_los=save_los,
        prng=prng,
    )


def project_photons_allsky(
    photon_prefix,
    event_prefix,
    normal,
    absorb_model=None,
    nH=None,
    abund_table="angr",
    no_shifting=False,
    center_vector=None,
    kernel="top_hat",
    save_los=False,
    prng=None,
):
    r"""
    Projects photons onto the sky sphere given a normal vector ("z" or "up" in
    spherical coordinates), and stores them in an HDF5 dataset which contains
    an event list.

    Parameters
    ----------
    photon_prefix : string
        The prefix of the filename(s) containing the photon list. If run in
        serial, the filename will be "{photon_prefix}.h5", if run in
        parallel, the filenames will be "{photon_prefix}.{mpi_rank}.h5".
    event_prefix : string
        The prefix of the filename(s) which will be written to contain the
        event list. If run in serial, the filename will be "{event_prefix}.h5",
        if run in parallel, the filename will be "{event_prefix}.{mpi_rank}.h5".
    normal : array-like
        The vector determining the "z" or "up" vector for the spherical coordinate
        system for the all-sky projection, something like [1.0, 2.0, -3.0]. It
        will be normalized before use.
    absorb_model : string
        A model for foreground galactic absorption, to simulate the
        absorption of events before being detected. Known options are "wabs"
        and "tbabs".
    nH : float, optional
        The foreground column density in units of 10^22 cm^{-2}. Only used
        if absorption is applied.
    abund_table : string
        The abundance table to be used for abundances in the
        absorption model (only used for TBabs). Default is set in the SOXS
        configuration file, the default for which is "angr".
        Built-in options are:
        "angr" : from Anders E. & Grevesse N. (1989, Geochimica et
        Cosmochimica Acta 53, 197)
        "aspl" : from Asplund M., Grevesse N., Sauval A.J. & Scott
        P. (2009, ARAA, 47, 481)
        "feld" : from Feldman U. (1992, Physica Scripta, 46, 202)
        "wilm" : from Wilms, Allen & McCray (2000, ApJ 542, 914
        except for elements not listed which are given zero abundance)
        "lodd" : from Lodders, K (2003, ApJ 591, 1220)
        "cl17.03" : the default abundance table in Cloudy 17.03
    no_shifting : boolean, optional
        If set, the photon energies will not be Doppler shifted. Default: False
    center_vector : a sequence of floats
        A vector defining what direction will be placed at the center of
        the lat/lon coordinate system. If not set, an arbitrary
        grid-aligned center_vector perpendicular to the normal is chosen.
    kernel : string, optional
        The kernel used when smoothing positions of X-rays originating from
        SPH particles, "gaussian" or "top_hat". Default: "top_hat".
    save_los : boolean, optional
        If True, save the line-of-sight radii in units of kpc to the events
        file. Default: False
    prng : integer or :class:`~numpy.random.RandomState` object
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers,
        such as for a test. Default is to use the :mod:`numpy.random`
        module.

    Returns
    -------
    A integer for the number of events created

    Examples
    --------
    >>> L = np.array([0.1,-0.2,0.3])
    >>> n_events = pyxsim.project_photons_allsky("my_photons.h5", "my_events.h5", L)
    """
    return _project_photons(
        "internal",
        photon_prefix,
        event_prefix,
        normal,
        [0.0, 0.0],
        absorb_model=absorb_model,
        nH=nH,
        abund_table=abund_table,
        no_shifting=no_shifting,
        kernel=kernel,
        save_los=save_los,
        north_vector=center_vector,
        prng=prng,
    )


class PhotonList:
    def __init__(self, filespec):
        """
        Read a PhotonList from disk to get information about it
        or to export to other formats.

        Parameters
        ----------
        filespec : str
            A filename, list of filenames, or globbable path
            that yields a single or list of HDF5 files containing
            event data.
        """
        from glob import glob

        if filespec.endswith(".h5"):
            filenames = glob(filespec)
        elif isinstance(filespec, list):
            if not np.all([fn.endswith(".h5") for fn in filespec]):
                raise RuntimeError("Not all filenames are valid!")
            filenames = filespec
        else:
            raise RuntimeError(f"Invalid PhotonList file spec: {filespec}")
        self.filenames = filenames
        self.filenames.sort()
        self.parameters = {}
        self.info = {}
        self.num_photons = []
        for i, fn in enumerate(self.filenames):
            with h5py.File(fn, "r") as f:
                p = f["parameters"]
                info = f["info"]
                self.num_photons.append(f["data"]["energy"].size)
                if i == 0:
                    for field in p:
                        if isinstance(p[field][()], (str, bytes)):
                            self.parameters[field] = p[field].asstr()[()]
                        else:
                            self.parameters[field] = p[field][()]
                    for k, v in info.attrs.items():
                        self.info[k] = v
        self.tot_num_photons = np.sum(self.num_photons)
        self.observer = self.parameters.get("observer", "external")

    def write_spectrum(self, specfile, emin, emax, nchan, overwrite=False):
        """
        Bin photon energies into a spectrum and write it to a FITS binary
        table. This is for an *unconvolved* spectrum.

        Parameters
        ----------
        specfile : string
            The name of the FITS file to be written.
        emin : float
            The minimum energy of the spectral bins in keV.
        emax : float
            The maximum energy of the spectral bins in keV.
        nchan : integer
            The number of channels.
        overwrite : boolean, optional
            Set to True to overwrite a previous file.
        """
        from astropy.io import fits

        spec = np.zeros(nchan)
        ebins = np.linspace(emin, emax, nchan + 1, endpoint=True)
        emid = 0.5 * (ebins[1:] + ebins[:-1])

        for fn in self.filenames:
            with h5py.File(fn, "r") as f:
                d = f["data"]
                spec += np.histogram(d["energy"][:], bins=ebins)[0]

        col1 = fits.Column(
            name="CHANNEL", format="1J", array=np.arange(nchan).astype("int32") + 1
        )
        col2 = fits.Column(name="ENERGY", format="1D", array=emid.astype("float64"))
        col3 = fits.Column(name="COUNTS", format="1J", array=spec.astype("int32"))
        col4 = fits.Column(
            name="COUNT_RATE", format="1D", array=spec / self.parameters["fid_exp_time"]
        )

        coldefs = fits.ColDefs([col1, col2, col3, col4])

        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = "SPECTRUM"

        tbhdu.header["DETCHANS"] = spec.shape[0]
        tbhdu.header["TOTCTS"] = spec.sum()
        tbhdu.header["EXPOSURE"] = self.parameters["fid_exp_time"]
        tbhdu.header["LIVETIME"] = self.parameters["fid_exp_time"]
        tbhdu.header["CONTENT"] = "pi"
        tbhdu.header["HDUCLASS"] = "OGIP"
        tbhdu.header["HDUCLAS1"] = "SPECTRUM"
        tbhdu.header["HDUCLAS2"] = "TOTAL"
        tbhdu.header["HDUCLAS3"] = "TYPE:I"
        tbhdu.header["HDUCLAS4"] = "COUNT"
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header["HDUVERS1"] = "1.1.0"
        tbhdu.header["CHANTYPE"] = "pi"
        tbhdu.header["BACKFILE"] = "none"
        tbhdu.header["CORRFILE"] = "none"
        tbhdu.header["POISSERR"] = True
        tbhdu.header["RESPFILE"] = "none"
        tbhdu.header["ANCRFILE"] = "none"
        tbhdu.header["MISSION"] = "none"
        tbhdu.header["TELESCOP"] = "none"
        tbhdu.header["INSTRUME"] = "none"
        tbhdu.header["AREASCAL"] = 1.0
        tbhdu.header["CORRSCAL"] = 0.0
        tbhdu.header["BACKSCAL"] = 1.0

        hdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])

        hdulist.writeto(specfile, overwrite=overwrite)
