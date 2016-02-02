import numpy as np
from yt.utilities.physical_ratios import \
    K_per_keV, mass_hydrogen_grams
from yt.frontends.stream.api import load_uniform_grid
from numpy.random import RandomState

my_prng = RandomState(24)

class BetaModelSource(object):
    def __init__(self):

        self.prng = RandomState(24)

        # Gas parameters
        R = 1.0 # Mpc
        r_c = 0.05 # Mpc
        rho_c = 0.04*mass_hydrogen_grams # g/cm**3
        beta = 1.
        self.kT = 6.0 # keV
        self.v_shift = 4.0e7 # cm/s
        self.v_width = 4.0e7 # cm/s

        # Dark matter parameters
        r_s = 0.350 # Mpc
        rho_s = 9.0e-26 # g/cm**3

        nx = 128
        ddims = (nx,nx,nx)

        x, y, z = np.mgrid[-R:R:nx*1j,
                           -R:R:nx*1j,
                           -R:R:nx*1j]

        r = np.sqrt(x**2+y**2+z**2)

        dens = np.zeros(ddims)
        dens[r <= R] = rho_c*(1.+(r[r <= R]/r_c)**2)**(-1.5*beta)
        dens[r > R] = 0.0
        pden = np.zeros(ddims)
        x = r[r <= R]/r_s
        pden[r <= R] = rho_s/(x*(1.+x)**2)
        pden[r > R] = 0.0
        temp = self.kT*K_per_keV*np.ones(ddims)
        bbox = np.array([[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
        velz = self.prng.normal(loc=self.v_shift,scale=self.v_width,size=ddims)
        dm_disp = 1000.*np.ones(ddims) # km/s

        data = {}
        data["density"] = (dens, "g/cm**3")
        data["dark_matter_density"] = (pden, "g/cm**3")
        data["dark_matter_dispersion"] = (dm_disp, "km/s")
        data["temperature"] = (temp, "K")
        data["velocity_x"] = (np.zeros(ddims), "cm/s")
        data["velocity_y"] = (np.zeros(ddims), "cm/s")
        data["velocity_z"] = (velz, "cm/s")

        self.ds = load_uniform_grid(data, ddims, length_unit=(2*R, "Mpc"),
                                    nprocs=64, bbox=bbox)

