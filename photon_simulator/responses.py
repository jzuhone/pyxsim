"""
Classes for response files. These generally will not be used
by the end-user, but will be employed by methods of EventList.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np
from yt.utilities.on_demand_imports import _astropy
from yt.funcs import mylog
from yt.units.yt_array import YTArray

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
        f = _astropy.pyfits.open(filename)
        self.elo = YTArray(f["SPECRESP"].data.field("ENERG_LO"), "keV")
        self.ehi = YTArray(f["SPECRESP"].data.field("ENERG_HI"), "keV")
        self.emid = 0.5*(self.elo+self.ehi)
        self.eff_area = YTArray(np.nan_to_num(f["SPECRESP"].data.field("SPECRESP")), "cm**2")
        self.filename = filename
        if rmffile is not None:
            rmf = RedistributionMatrixFile(rmffile)
            self.eff_area *= rmf.weights
        else:
            mylog.warning("You specified an ARF but not an RMF. This is ok if you know "
                          "a priori that the responses are normalized properly. If not, "
                          "you may get inconsistent results.")
        f.close()

    def __str__(self):
        return self.filename

    def detect_events(self, energy, prng=np.random):
        """
        Use the ARF to determine a subset of photons which will be
        detected. Returns a boolean NumPy array which is the same
        is the same size as the number of photons, wherever it is 
        "true" means those photons have been detected.

        Parameters
        ----------
        energy : YTArray
            The energies of the photons to attempt to detect.
        prng : NumPy `RandomState` object or numpy.random
            A pseudo-random number generator. Typically will only be specified if you
            have a reason to generate the same set of random numbers, such as for a
            test. Default is the numpy.random module.
        """
        earea = np.interp(energy, self.emid, self.eff_area, left=0.0, right=0.0)
        randvec = self.max_area*prng.uniform(size=energy.shape)
        return randvec < earea

    @property
    def max_area(self):
        return self.eff_area.max()

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
        mylog.info("Reading response matrix file (RMF): %s" % filename)
        self.handle = _astropy.pyfits.open(filename)
        if "MATRIX" in self.handle:
            self.mat_key = "MATRIX"
        elif "SPECRESP MATRIX" in self.handle:
            self.mat_key = "SPECRESP MATRIX"
        else:
            raise RuntimeError("Cannot find the response matrix in the RMF "
                               "file %s! " % filename+"It should be named "
                                                      "\"MATRIX\" or \"SPECRESP MATRIX\".")
        self.data = self.handle[self.mat_key].data
        self.header = self.handle[self.mat_key].header
        self.ebounds = self.handle["EBOUNDS"].data
        self.weights = np.array([w.sum() for w in self.data["MATRIX"]])
        self.filename = filename

    def __str__(self):
        return self.filename
