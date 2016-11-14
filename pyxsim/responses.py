"""
Classes for response files. These generally will not be used
by the end-user, but will be employed by methods of InstrumentSimulators.
"""

import numpy as np
from yt.utilities.on_demand_imports import _astropy
from yt.units.yt_array import YTArray
from pyxsim.utils import mylog, check_file_location

class AuxiliaryResponseFile(object):
    r"""
    A class for auxiliary response files (ARFs).

    Parameters
    ----------
    filename : string
        The filename of the ARF to be read.
    rmffile : string
        The filename of a redistribution matrix
        file (RMF) to be read. This is sometimes
        needed because some ARFs are unnormalized.

    Examples
    --------
    >>> arf = AuxiliaryResponseFile("sxt-s_120210_ts02um_intallpxl.arf",
    ...                             rmffile="ah_sxs_5ev_basefilt_20100712.rmf")
    """
    def __init__(self, filename, rmffile=None):
        self.filename = check_file_location(filename, "response_files")
        f = _astropy.pyfits.open(self.filename)
        self.elo = YTArray(f["SPECRESP"].data.field("ENERG_LO"), "keV")
        self.ehi = YTArray(f["SPECRESP"].data.field("ENERG_HI"), "keV")
        self.emid = 0.5*(self.elo+self.ehi)
        self.eff_area = YTArray(np.nan_to_num(f["SPECRESP"].data.field("SPECRESP")), "cm**2")
        if rmffile is not None:
            rmf = RedistributionMatrixFile(rmffile)
            self.eff_area *= rmf.weights
            self.rmffile = rmf.filename
        else:
            mylog.warning("You specified an ARF but not an RMF. This is ok if you know "
                          "a priori that the responses are normalized properly. If not, "
                          "you may get inconsistent results.")
        f.close()

    def __str__(self):
        return self.filename

    def detect_events(self, energy, area, prng=None):
        """
        Use the ARF to determine a subset of photons which will be
        detected. Returns a boolean NumPy array which is the same
        is the same size as the number of photons, wherever it is 
        "true" means those photons have been detected.

        Parameters
        ----------
        energy : :class:`~yt.units.yt_array.YTArray`
            The energies of the photons to attempt to detect, in keV.
        area : float, tuple, or YTQuantity
            The collecting area associated with the event energies. If a floating-point
            number, is assumed to be in cm^2. 
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`~numpy.random` module.
        """
        if prng is None:
            prng = np.random
        earea = np.interp(energy, self.emid, self.eff_area, left=0.0, right=0.0)
        randvec = area.v*prng.uniform(size=energy.shape)
        return randvec < earea

    @property
    def max_area(self):
        return self.eff_area.max()

    def interpolate_area(self, energy):
        """
        Interpolate the effective area to the energies provided by the supplied *energy* array.
        """
        earea = np.interp(energy, self.emid, self.eff_area, left=0.0, right=0.0)
        return YTArray(earea, "cm**2")

class RedistributionMatrixFile(object):
    r"""
    A class for redistribution matrix files (RMFs).

    Parameters
    ----------
    filename : string
        The filename of the RMF to be read.

    Examples
    --------
    >>> rmf = RedistributionMatrixFile("acisi_aimpt_cy17.rmf")
    """
    def __init__(self, filename):
        self.filename = check_file_location(filename, "response_files")
        self.handle = _astropy.pyfits.open(self.filename)
        names = [h.name.strip() for h in self.handle]
        if "MATRIX" in names:
            self.mat_key = "MATRIX"
        elif "SPECRESP MATRIX" in names:
            self.mat_key = "SPECRESP MATRIX"
        else:
            raise RuntimeError("Cannot find the response matrix in the RMF "
                               "file %s! " % filename+"It should be named "
                                                      "\"MATRIX\" or \"SPECRESP MATRIX\".")
        self.data = self.handle[self.mat_key].data
        self.header = self.handle[self.mat_key].header
        self.num_mat_columns = len(self.handle[self.mat_key].columns)
        self.ebounds = self.handle["EBOUNDS"].data
        self.ebounds_header = self.handle["EBOUNDS"].header
        self.weights = np.array([w.sum() for w in self.data["MATRIX"]])
        self.elo = self.data["ENERG_LO"]
        self.ehi = self.data["ENERG_HI"]
        self.n_de = self.elo.size
        self.n_ch = self.ebounds["CHANNEL"].size

        num = 0
        for i in range(1, self.num_mat_columns+1):
            if self.header["TTYPE%d" % i] == "F_CHAN":
                num = i
                break
        self.cmin = self.header["TLMIN%d" % num]
        self.cmax = self.header["TLMAX%d" % num]

    def __str__(self):
        return self.filename
