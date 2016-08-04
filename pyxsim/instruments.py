import numpy as np
from pyxsim.event_list import EventList
from pyxsim.responses import AuxiliaryResponseFile, \
    RedistributionMatrixFile
from pyxsim.utils import mylog
from yt.funcs import get_pbar, ensure_numpy_array, \
    iterable
from yt.units.yt_array import YTQuantity, YTArray
from yt.utilities.on_demand_imports import _astropy
from copy import deepcopy

sigma_to_fwhm = 2.*np.sqrt(2.*np.log(2.))

class InstrumentSimulator(object):
    def __init__(self, dtheta, nx, psf_scale, arf,
                 rmf):
        """
        Construct an instrument simulator.

        Parameters
        ----------
        dtheta : float
            The width of the central pixel in degrees. 
        nx : integer
            The number of resolution elements on a side across
            the field of view.
        psf_scale : float
            The FWHM of the Gaussian PSF in degrees.
        arf : string
            The path to the ARF file that will be used for
            the effective area.
        rmf : string
            The path to the RMF file that will be used for
            the spectral response matrix. 

        Examples
        --------
        >>> from pyxsim import InstrumentSimulator
        >>> ACIS_S = InstrumentSimulator(0.0001366667, 8192, 0.0001388889,
        ...                              "aciss_aimpt_cy18.arf",
        ...                              "aciss_aimpt_cy18.rmf")
        """
        self.dtheta = dtheta
        self.nx = nx
        self.psf_scale = psf_scale
        self.arf = arf
        self.rmf = rmf

    def __call__(self, events, rebin=True,
                 convolve_psf=True, convolve_arf=True, 
                 convolve_rmf=True, prng=None):
        new_events = EventList(deepcopy(events.events), 
                               events.parameters.copy(), events.wcs.copy())
        if prng is None:
            prng = np.random
        if rebin:
            self.rebin(new_events)
        if convolve_psf:
            self.convolve_with_psf(new_events, prng)
        if convolve_arf:
            new_events["xsky"]
            new_events["ysky"]
            self.apply_effective_area(new_events, prng)
            if convolve_rmf:
                self.convolve_energies(new_events, prng)
        return new_events

    def rebin(self, events):
        """
        Rebin event positions to a new binning with the same celestial
        coordinates.
        """
        new_wcs = _astropy.pywcs.WCS(naxis=2)
        new_wcs.wcs.crval = events.parameters["sky_center"].d
        new_wcs.wcs.crpix = np.array([0.5*(self.nx+1)]*2)
        new_wcs.wcs.cdelt = [-self.dtheta, self.dtheta]
        new_wcs.wcs.ctype = ["RA---TAN","DEC--TAN"]
        new_wcs.wcs.cunit = ["deg"]*2
        xpix, ypix = new_wcs.wcs_world2pix(events["xsky"], events["ysky"], 1)
        events.events['xpix'] = xpix
        events.events['ypix'] = ypix
        xsky, ysky = new_wcs.wcs_pix2world(events["xpix"], events["ypix"], 1)
        events.events['xsky'] = xsky
        events.events['ysky'] = ysky
        events.parameters['pix_center'] = new_wcs.wcs.crpix[:]
        events.parameters['dtheta'] = YTQuantity(self.dtheta, "deg")
        events.wcs = new_wcs

    def convolve_with_psf(self, events, prng):
        r"""
        Convolve the event positions with a PSF.
        """
        dtheta = events.parameters["dtheta"]
        psf = lambda n: prng.normal(scale=self.psf_scale/sigma_to_fwhm/dtheta, size=n)
        events.events["xpix"] += psf(events.num_events)
        events.events["ypix"] += psf(events.num_events)
        xsky, ysky = events.wcs.wcs_pix2world(events["xpix"], events["ypix"], 1)
        events.events['xsky'] = xsky
        events.events['ysky'] = ysky

    def apply_effective_area(self, events, prng):
        """
        Convolve the events with a ARF file.
        """
        mylog.info("Applying energy-dependent effective area.")
        arf = AuxiliaryResponseFile(self.arf, rmffile=self.rmf)
        # If the area which was used to create the events is smaller than
        # the maximum area in the ARF, scream loudly.
        if events.parameters["Area"] < arf.max_area:
            raise RuntimeError("The area used to create the events is less than "
                               "the maximum of the effective area curve! Re-create the "
                               "events with a collecting area higher than %s!" % arf.max_area)
        detected = arf.detect_events(events["eobs"], events.parameters["Area"], prng=prng)
        mylog.info("%s events detected." % detected.sum())
        for key in ["xpix", "ypix", "xsky", "ysky", "eobs"]:
            events.events[key] = events[key][detected]
        events.parameters["ARF"] = arf.filename
        events.num_events = len(events.events["eobs"])

    def convolve_energies(self, events, prng):
        """
        Convolve the events with a RMF file.
        """
        mylog.info("Reading response matrix file (RMF): %s" % self.rmf)
        rmf = RedistributionMatrixFile(self.rmf)

        elo = rmf.data["ENERG_LO"]
        ehi = rmf.data["ENERG_HI"]
        n_de = elo.shape[0]
        mylog.info("Number of energy bins in RMF: %d" % n_de)
        mylog.info("Energy limits: %g %g" % (min(elo), max(ehi)))

        n_ch = len(rmf.ebounds["CHANNEL"])
        mylog.info("Number of channels in RMF: %d" % n_ch)

        eidxs = np.argsort(events["eobs"])
        sorted_e = events["eobs"][eidxs].d

        detectedChannels = []

        # run through all photon energies and find which bin they go in
        fcurr = 0
        last = sorted_e.shape[0]

        pbar = get_pbar("Scattering energies with RMF", last)

        for (k, low), high in zip(enumerate(elo), ehi):
            # weight function for probabilities from RMF
            weights = np.nan_to_num(np.float64(rmf.data["MATRIX"][k]))
            weights /= weights.sum()
            # build channel number list associated to array value,
            # there are groups of channels in rmfs with nonzero probabilities
            trueChannel = []
            f_chan = ensure_numpy_array(np.nan_to_num(rmf.data["F_CHAN"][k]))
            n_chan = ensure_numpy_array(np.nan_to_num(rmf.data["N_CHAN"][k]))
            if not iterable(f_chan):
                f_chan = [f_chan]
                n_chan = [n_chan]
            for start, nchan in zip(f_chan, n_chan):
                if nchan == 0:
                    trueChannel.append(start)
                else:
                    trueChannel += list(range(start, start+nchan))
            if len(trueChannel) > 0:
                for q in range(fcurr, last):
                    if low <= sorted_e[q] < high:
                        channelInd = prng.choice(len(weights), p=weights)
                        fcurr += 1
                        pbar.update(fcurr)
                        detectedChannels.append(trueChannel[channelInd])
                    else:
                        break
        pbar.finish()

        for key in ["xpix", "ypix", "xsky", "ysky"]:
            events.events[key] = events[key][eidxs]

        events.events["eobs"] = YTArray(sorted_e, "keV")
        events.events[rmf.header["CHANTYPE"]] = np.array(detectedChannels, dtype="int")

        events.parameters["RMF"] = rmf.filename
        events.parameters["ChannelType"] = rmf.header["CHANTYPE"]
        events.parameters["Telescope"] = rmf.header["TELESCOP"]
        events.parameters["Instrument"] = rmf.header["INSTRUME"]
        events.parameters["Mission"] = rmf.header.get("MISSION","")

# Specific instrument approximations

ACIS_S = InstrumentSimulator(0.0001366667, 8192, 0.0001388889,
                             "aciss_aimpt_cy18.arf",
                             "aciss_aimpt_cy18.rmf")
ACIS_I = InstrumentSimulator(0.0001366667, 8192, 0.0001388889,
                             "acisi_aimpt_cy18.arf",
                             "acisi_aimpt_cy18.rmf")
Hitomi_SXS = InstrumentSimulator(8.512516E-03, 6, 0.02,
                                 "sxt-s_140505_ts02um_intallpxl.arf",
                                 "ah_sxs_5ev_20130806.rmf")
Athena_WFI = InstrumentSimulator(6.207043E-04, 1024, 0.001388888888888889,
                                 "athena_wfi_1469_onaxis_w_filter_v20150326.arf", 
                                 "athena_wfi_rmf_v20150326.rmf")
Athena_XIFU = InstrumentSimulator(1.265282E-03, 66, 0.001388888888888889,
                                  "athena_xifu_1469_onaxis_pitch265um_v20150327.arf",
                                  "athena_xifu_rmf_v20150327.rmf")
XRS_Imager = InstrumentSimulator(9.167325E-05, 4096, 0.0001388889,
                                 "xrs_hdxi.arf", "xrs_hdxi.rmf")
XRS_Calorimeter = InstrumentSimulator(0.0002864789, 300, 0.0001388889,
                                      "xrs_calorimeter.arf",
                                      "xrs_calorimeter.rmf")
