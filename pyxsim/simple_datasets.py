from yt.frontends.stream.api import load_amr_grids, \
    load_uniform_grid, load_particles, refine_amr
import numpy as np

def create_spherical_dataset(profiles, ddims, width, nprocs=64):
    R = 0.5*width
    nx, ny, nz = ddims
    x, y, z = np.mgrid[-R:R:nx * 1j,
                       -R:R:ny * 1j,
                       -R:R:nz * 1j]

    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    dens = np.zeros(ddims)
    dens[r <= R] = rho_c * (1. + (r[r <= R] / r_c) ** 2) ** (-1.5 * beta)
    dens[r > R] = 0.0
    pden = np.zeros(ddims)
    x = r[r <= R] / r_s
    pden[r <= R] = rho_s / (x * (1. + x) ** 2)
    pden[r > R] = 0.0
    temp = self.kT * K_per_keV * np.ones(ddims)
    bbox = np.array([[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]])
    velz = self.prng.normal(loc=v_shift, scale=v_width, size=ddims)
    dm_disp = 1000. * np.ones(ddims)  # km/s

    data = {}
    data["density"] = (dens, "g/cm**3")
    data["dark_matter_density"] = (pden, "g/cm**3")
    data["dark_matter_dispersion"] = (dm_disp, "km/s")
    data["temperature"] = (temp, "K")
    data["velocity_x"] = (np.zeros(ddims), "cm/s")
    data["velocity_y"] = (np.zeros(ddims), "cm/s")
    data["velocity_z"] = (velz, "cm/s")

    ds = load_uniform_grid(data, ddims, nprocs=nprocs)

    return ds

def create_spherical_amr_dataset(profiles, ddims, ref_level, width, nprocs=64):
    data = {}
    ds = load_amr_grids()
    return ds

def create_spherical_particle_dataset():
    pass