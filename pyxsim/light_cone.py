import time

from soxs.utils import parse_prng
from yt.loaders import load
from yt.utilities.parallel_tools.parallel_analysis_interface import communication_system
from yt_astro_analysis.cosmological_observation.api import LightCone

from pyxsim.photon_list import make_photons, project_photons
from pyxsim.utils import parse_value

comm = communication_system.communicators[-1]

axes_lookup = [(1, 2), (2, 0), (0, 1)]


class XrayLightCone(LightCone):
    def __init__(
        self,
        parameter_filename,
        simulation_type,
        near_redshift,
        far_redshift,
        seed=None,
        use_minimum_datasets=True,
        deltaz_min=0.0,
        minimum_coherent_box_fraction=0.0,
    ):
        if seed is None:
            seed = time.time()
        super(XrayLightCone, self).__init__(
            parameter_filename,
            simulation_type,
            near_redshift,
            far_redshift,
            use_minimum_datasets=use_minimum_datasets,
            deltaz_min=deltaz_min,
            minimum_coherent_box_fraction=minimum_coherent_box_fraction,
        )
        self.calculate_light_cone_solution(seed=seed)

    def generate_events(
        self,
        photon_prefix,
        event_prefix,
        area,
        exp_time,
        angular_width,
        source_model,
        sky_center,
        parameters=None,
        velocity_fields=None,
        absorb_model=None,
        nH=None,
        no_shifting=False,
        sigma_pos=None,
        prng=None,
    ):
        """
        Generate projected events from a light cone simulation.

        Parameters
        ----------
        photon_prefix : string
            The prefix of the filename(s) containing the photon list. If run in
            serial, the filename will be "{photon_prefix}.lc{i}.h5", where i
            iterates over the elements of the light cone solution. If run in
            parallel, the filenames will be "{photon_prefix}.lc{i}.{mpi_rank}.h5".
        event_prefix : string
            The prefix of the filename(s) which will be written to contain the
            event list. If run in serial, the filename will be "{event_prefix}.h5",
            if run in parallel, the filename will be "{event_prefix}.{mpi_rank}.h5".
        area : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
            The collecting area to determine the number of events. If units are
            not specified, it is assumed to be in cm^2.
        exp_time : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
            The exposure time to determine the number of events. If units are
            not specified, it is assumed to be in seconds.
        angular_width : float, (value, unit) tuple, or :class:`~yt.units.yt_array.YTQuantity`
            The angular width of the light cone simulation. If units are not
            specified, it is assumed to be in degrees.
        source_model : :class:`~pyxsim.source_models.sources.SourceModel`
            A source model used to generate the events.
        sky_center : array-like
            Center RA, Dec of the events in degrees.
        parameters : dict, optional
            A dictionary of parameters to be passed for the source model to use,
            if necessary.
        velocity_fields : list of fields
            The yt fields to use for the velocity. If not specified, the following will
            be assumed:
            ['velocity_x', 'velocity_y', 'velocity_z'] for grid datasets
            ['particle_velocity_x', 'particle_velocity_y', 'particle_velocity_z'] for particle datasets
        absorb_model : string
            A model for foreground galactic absorption, to simulate the absorption
            of events before being detected. Known options are "wabs" and "tbabs".
        nH : float, optional
            The foreground column density in units of 10^22 cm^{-2}. Only used if
            absorption is applied.
        no_shifting : boolean, optional
            If set, the photon energies will not be Doppler shifted.
        sigma_pos : float, optional
            Apply a gaussian smoothing operation to the sky positions of the
            events. This may be useful when the binned events appear blocky due
            to their uniform distribution within simulation cells. However, this
            will move the events away from their originating position on the
            sky, and so may distort surface brightness profiles and/or spectra.
            Should probably only be used for visualization purposes. Supply a
            float here to smooth with a standard deviation with this fraction
            of the cell size. Default: None
        prng : integer or :class:`~numpy.random.RandomState` object
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is to use the :mod:`numpy.random` module.
        """
        prng = parse_prng(prng)

        area = parse_value(area, "cm**2")
        exp_time = parse_value(exp_time, "s")
        aw = parse_value(angular_width, "deg")

        for i, output in enumerate(self.light_cone_solution):
            ds = load(output["filename"])
            dw = ds.domain_width
            ax = output["projection_axis"]
            c = output["projection_center"] * dw + ds.domain_left_edge
            le = c.copy()
            re = c.copy()
            width = ds.quan(aw * output["box_width_per_angle"], "unitary").to(
                "code_length"
            )
            depth = dw[ax].to("code_length") * output["box_depth_fraction"]
            le[ax] -= 0.5 * depth
            re[ax] += 0.5 * depth
            for off_ax in axes_lookup[ax]:
                le[off_ax] -= 0.5 * width
                re[off_ax] += 0.5 * width
            reg = ds.box(le, re)
            pprefix = f"{photon_prefix}.lc{i}"
            n_photons, n_cells = make_photons(
                pprefix,
                reg,
                output["redshift"],
                area,
                exp_time,
                source_model,
                parameters=parameters,
                center=c,
                velocity_fields=velocity_fields,
                cosmology=ds.cosmology,
            )
            eprefix = f"{event_prefix}.lc{i}"
            project_photons(
                pprefix,
                eprefix,
                "xyz"[ax],
                sky_center,
                absorb_model=absorb_model,
                nH=nH,
                no_shifting=no_shifting,
                sigma_pos=sigma_pos,
                prng=prng,
            )

            comm.barrier()
