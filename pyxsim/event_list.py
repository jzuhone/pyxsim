"""
Classes for generating lists of detected events
"""
from six import string_types
from collections import defaultdict
import numpy as np
from yt.funcs import mylog, get_pbar, ensure_numpy_array, \
    iterable, ensure_list
from yt.utilities.fits_image import assert_same_wcs
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.units.yt_array import YTQuantity, YTArray, uconcatenate
from yt.utilities.on_demand_imports import _astropy
import warnings
import os
import h5py
from pyxsim.utils import force_unicode, validate_parameters
from pyxsim.responses import AuxiliaryResponseFile

class EventList(object):

    def __init__(self, events, parameters):
        self.events = events
        self.parameters = parameters
        self.num_events = events["xpix"].shape[0]
        self.wcs = _astropy.pywcs.WCS(naxis=2)
        self.wcs.wcs.crpix = parameters["pix_center"]
        self.wcs.wcs.crval = parameters["sky_center"].d
        self.wcs.wcs.cdelt = [-parameters["dtheta"].value, parameters["dtheta"].value]
        self.wcs.wcs.ctype = ["RA---TAN","DEC--TAN"]
        self.wcs.wcs.cunit = ["deg"]*2

    def keys(self):
        return self.events.keys()

    def has_key(self, key):
        return key in self.keys()

    def items(self):
        return self.events.items()

    def values(self):
        return self.events.values()

    def __getitem__(self,key):
        if key not in self.events:
            if key == "xsky" or key == "ysky":
                x,y = self.wcs.wcs_pix2world(self.events["xpix"], self.events["ypix"], 1)
                self.events["xsky"] = YTArray(x, "degree")
                self.events["ysky"] = YTArray(y, "degree")
        return self.events[key]

    def __repr__(self):
        return self.events.__repr__()

    def __contains__(self, key):
        return key in self.events

    def __add__(self, other):
        assert_same_wcs(self.wcs, other.wcs)
        validate_parameters(self.parameters, other.parameters)
        events = {}
        for item1, item2 in zip(self.items(), other.items()):
            k1, v1 = item1
            k2, v2 = item2
            events[k1] = uconcatenate([v1,v2])
        return EventList(events, self.parameters)

    def convolve_energies(self, rmf, prng=np.random, clobber_channels=False):
        """
        Convolve the events with a RMF file.

        Parameters
        ----------
        rmf : :class:`~pyxsim.responses.RedistributionMatrixFile`
            The RMF to use when convolving the energies.
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        clobber_channels : boolean, optional
            If channels have already been determined, set this to True to
            re-calculate the channels.

        Examples
        --------
        >>> from numpy.random import RandomState
        >>> prng = RandomState(25)
        >>> rmf = RedistributionMatrixFile("pn-med.rmf")
        >>> events.convolve_energies(rmf, prng=prng, clobber_channels=True)
        """
        if "RMF" in self.parameters and rmf.filename != self.parameters["RMF"]:
            err = "This EventList is already associated with an RMF: %s," % self.parameters["RMF"]
            err += "but you want to convolve with a different RMF: %s!" % rmf.filename
            raise RuntimeError(err)

        if ("PI" in self or "PHA" in self) and not clobber_channels:
            raise RuntimeError("You've already convolved these events with "
                               "an RMF! If you want to overwrite them, set "
                               "clobber_channels=True!")

        mylog.info("Reading response matrix file (RMF): %s" % rmf.filename)

        elo = rmf.data["ENERG_LO"]
        ehi = rmf.data["ENERG_HI"]
        n_de = elo.shape[0]
        mylog.info("Number of energy bins in RMF: %d" % n_de)
        mylog.info("Energy limits: %g %g" % (min(elo), max(ehi)))

        n_ch = len(rmf.ebounds["CHANNEL"])
        mylog.info("Number of channels in RMF: %d" % n_ch)

        eidxs = np.argsort(self.events["eobs"])
        sorted_e = self.events["eobs"][eidxs].d

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


        self.events["xpix"] = self.events["xpix"][eidxs]
        self.events["ypix"] = self.events["ypix"][eidxs]
        if "xsky" in self.events:
            self.events["xsky"] = self.events["xsky"][eidxs]
        if "ysky" in self.events:
            self.events["ysky"] = self.events["ysky"][eidxs]
        self.events["eobs"] = YTArray(sorted_e, "keV")
        self.events[rmf.header["CHANTYPE"]] = np.array(detectedChannels, dtype="int")

        self.parameters["RMF"] = rmf.filename
        self.parameters["ChannelType"] = rmf.header["CHANTYPE"]
        self.parameters["Telescope"] = rmf.header["TELESCOP"]
        self.parameters["Instrument"] = rmf.header["INSTRUME"]
        self.parameters["Mission"] = rmf.header.get("MISSION","")

    def add_point_sources(self, positions, energy_bins, spectra,
                          prng=None, absorb_model=None):
        r"""
        Add point source events to an :class:`~photon_simulator.event_list.EventList`.

        Parameters
        ----------
        positions : array of pixel positions, shape 2xN
            The positions of the point sources in pixel space, where N is the
            number of point sources
        energy_bins : :class:`~yt.units.yt_array.YTArray` with units of keV, shape M+1
            The edges of the energy bins for the spectra, where M is the number of
            bins
        spectra : list (size N) of :class:`~yt.units.yt_array.YTArray`s with units of photons/s/cm^2, each with shape M
            The spectra for the point sources, where M is the number of bins and N is
            the number of point sources
        prng : :class:`~numpy.random.RandomState` object or :mod:`numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        absorb_model : :class:`~photon_simulator.spectral_models.TableAbsorbModel` or :class:`~photon_simulator.spectral_models.XSpecAbsorbModel`, optional
            A model for galactic absorption.
        """
        if prng is None:
            prng = np.random

        spectra = ensure_list(spectra)
        if positions.shape == (2,):
            positions = positions.reshape(2,1)

        x = [self.events["xpix"]]
        y = [self.events["ypix"]]
        e = [self.events["eobs"]]

        for pos, spectrum in zip(positions, spectra):
            eobs = self._add_events(energy_bins, spectrum, prng, absorb_model)
            x, y = self.wcs.wcs_world2pix(pos[0], pos[1], 1)
            ne = len(eobs)
            x.append([x]*ne)
            y.append([y]*ne)
            e.append(eobs)

        events = {}
        events["xpix"] = uconcatenate(x)
        events["ypix"] = uconcatenate(y)
        events["eobs"] = uconcatenate(e)

        return EventList(events, self.parameters)

    def add_background(self, energy_bins, spectrum,
                       prng=np.random, absorb_model=None):
        r"""
        Add background events to an :class:`~photon_simulator.event_list.EventList`.

        Parameters
        ----------
        energy_bins : :class:`~yt.units.yt_array.YTArray` with units of keV, size M+1
            The edges of the energy bins for the spectra, where M is the number of
            bins
        spectrum : :class:`~yt.units.yt_array.YTArray` with units of photons/s/cm^2, size M
            The spectrum for the background, where M is the number of bins.
        prng : :class:`~numpy.random.RandomState` object or :mod:`numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        absorb_model : :class:`~photon_simulator.spectral_models.TableAbsorbModel` or :class:`~photon_simulator.spectral_models.XSpecAbsorbModel`, optional
            A model for galactic absorption.
        """
        if prng is None:
            prng = np.random

        eobs = self._add_events(energy_bins, spectrum, prng, absorb_model)
        ne = len(eobs)
        x = np.random.uniform(low=0.5, high=2.*self.parameters["pix_center"][0]-0.5, size=ne)
        y = np.random.uniform(low=0.5, high=2.*self.parameters["pix_center"][1]-0.5, size=ne)

        events = {}
        events["xpix"] = uconcatenate([x, self.events["xpix"]])
        events["ypix"] = uconcatenate([y, self.events["ypix"]])
        events["eobs"] = uconcatenate([eobs, self.events["eobs"]])

        return EventList(events, self.parameters)

    def _add_events(self, ebins, spectrum, prng, absorb_model):
        exp_time = self.parameters["exp_time"]
        area = self.parameters["area"]
        flux = spectrum.sum()
        num_photons = np.uint64(exp_time*area*flux)
        cumspec = np.cumsum(spectrum)
        cumspec = np.insert(cumspec, 0, 0.0)
        cumspec /= cumspec[-1]
        randvec = prng.uniform(size=num_photons)
        randvec.sort()
        e = YTArray(np.interp(randvec, cumspec, ebins), "keV")

        if absorb_model is None:
            not_abs = np.ones(e.shape, dtype='bool')
        else:
            mylog.info("Absorbing.")
            absorb_model.prepare_spectrum()
            emid = absorb_model.emid
            aspec = absorb_model.get_spectrum()
            absorb = np.interp(e, emid, aspec, left=0.0, right=0.0)
            randvec = aspec.max()*prng.uniform(size=e.shape)
            not_abs = randvec < absorb
            absorb_model.cleanup_spectrum()

        if "ARF" in self.parameters:
            rmffile = self.parameters.get("RMF", None)
            arf = AuxiliaryResponseFile(self.parameters["ARF"],
                                        rmffile=rmffile)
            detected = arf.detect_events(e, prng=prng)
        else:
            detected = np.ones(e.shape, dtype='bool')

        detected = np.logical_and(not_abs, detected)

        return e[detected]

    def filter_events(self, region):
        """
        Filter events using a ds9 *region*. Requires the ``pyregion`` package.
        Returns a new :class:`~pyxsim.event_list.EventList`.
        """
        import pyregion
        import os
        if os.path.exists(region):
            reg = pyregion.open(region)
        else:
            reg = pyregion.parse(region)
        r = reg.as_imagecoord(header=self.wcs.to_header())
        f = r.get_filter()
        idxs = f.inside_x_y(self["xpix"], self["ypix"])
        if idxs.sum() == 0:
            raise RuntimeError("No events are inside this region!")
        new_events = {}
        for k, v in self.items():
            new_events[k] = v[idxs]
        return EventList(new_events, self.parameters)

    def convolve_with_psf(self, psf, prng=None):
        r"""
        Convolve the event positions with a PSF.

        Parameters
        ----------
        psf : float or function
            The PSF to convolve the photon positions with. If a function, it must take
            the number of events as an argument. If a float, it will be assumed to be the
            standard deviation of a Gaussian PSF.
        prng : :class:`~numpy.random.RandomState` object or :mod:`numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.

        Examples
        --------
        """
        if prng is None:
            prng = np.random
        dtheta = self.parameters['dtheta']
        new_events = self.events.deepcopy()
        if isinstance(psf, float):
            psf = lambda n: prng.normal(scale=psf/dtheta, size=n)
        new_events["xpix"] += psf(self.num_events)
        new_events["ypix"] += psf(self.num_events)
        return EventList(new_events, self.parameters)

    @classmethod
    def from_h5_file(cls, h5file):
        """
        Initialize an :class:`~pyxsim.event_list.EventList` from a HDF5 file with filename *h5file*.
        """
        events = {}
        parameters = {}

        f = h5py.File(h5file, "r")

        p = f["/parameters"]
        parameters["ExposureTime"] = YTQuantity(p["exp_time"].value, "s")
        parameters["Area"] = YTQuantity(p["area"].value, "cm**2")
        parameters["Redshift"] = p["redshift"].value
        parameters["AngularDiameterDistance"] = YTQuantity(p["d_a"].value, "Mpc")
        parameters["sky_center"] = YTArray(p["sky_center"][:], "deg")
        parameters["dtheta"] = YTQuantity(p["dtheta"].value, "deg")
        parameters["pix_center"] = p["pix_center"][:]
        if "rmf" in p:
            parameters["RMF"] = force_unicode(p["rmf"].value)
        if "arf" in p:
            parameters["ARF"] = force_unicode(p["arf"].value)
        if "channel_type" in p:
            parameters["ChannelType"] = force_unicode(p["channel_type"].value)
        if "mission" in p:
            parameters["Mission"] = force_unicode(p["mission"].value)
        if "telescope" in p:
            parameters["Telescope"] = force_unicode(p["telescope"].value)
        if "instrument" in p:
            parameters["Instrument"] = force_unicode(p["instrument"].value)

        d = f["/data"]
        events["xpix"] = d["xpix"][:]
        events["ypix"] = d["ypix"][:]
        events["eobs"] = YTArray(d["eobs"][:], "keV")
        if "pi" in d:
            events["PI"] = d["pi"][:]
        if "pha" in d:
            events["PHA"] = d["pha"][:]

        f.close()

        return cls(events, parameters)

    @classmethod
    def from_fits_file(cls, fitsfile):
        """
        Initialize an :class:`~pyxsim.event_list.EventList` from a FITS file with filename *fitsfile*.
        """
        hdulist = _astropy.pyfits.open(fitsfile)

        tblhdu = hdulist["EVENTS"]

        events = {}
        parameters = {}

        parameters["ExposureTime"] = YTQuantity(tblhdu.header["EXPOSURE"], "s")
        parameters["Area"] = YTQuantity(tblhdu.header["AREA"], "cm**2")
        parameters["Redshift"] = tblhdu.header["REDSHIFT"]
        parameters["AngularDiameterDistance"] = YTQuantity(tblhdu.header["D_A"], "Mpc")
        if "RMF" in tblhdu.header:
            parameters["RMF"] = tblhdu["RMF"]
        if "ARF" in tblhdu.header:
            parameters["ARF"] = tblhdu["ARF"]
        if "CHANTYPE" in tblhdu.header:
            parameters["ChannelType"] = tblhdu["CHANTYPE"]
        if "MISSION" in tblhdu.header:
            parameters["Mission"] = tblhdu["MISSION"]
        if "TELESCOP" in tblhdu.header:
            parameters["Telescope"] = tblhdu["TELESCOP"]
        if "INSTRUME" in tblhdu.header:
            parameters["Instrument"] = tblhdu["INSTRUME"]
        parameters["sky_center"] = YTArray([tblhdu["TCRVL2"],tblhdu["TCRVL3"]], "deg")
        parameters["pix_center"] = np.array([tblhdu["TCRVL2"],tblhdu["TCRVL3"]])
        parameters["dtheta"] = YTQuantity(tblhdu["TCRVL3"], "deg")
        events["xpix"] = tblhdu.data.field("X")
        events["ypix"] = tblhdu.data.field("Y")
        events["eobs"] = YTArray(tblhdu.data.field("ENERGY")/1000., "keV")
        if "PI" in tblhdu.columns.names:
            events["PI"] = tblhdu.data.field("PI")
        if "PHA" in tblhdu.columns.names:
            events["PHA"] = tblhdu.data.field("PHA")

        return cls(events, parameters)

    @parallel_root_only
    def write_fits_file(self, fitsfile, clobber=False):
        """
        Write events to a FITS binary table file with filename *fitsfile*.
        Set *clobber* to True if you need to overwrite a previous file.
        """
        from astropy.time import Time, TimeDelta
        pyfits = _astropy.pyfits

        exp_time = float(self.parameters["ExposureTime"])

        t_begin = Time.now()
        dt = TimeDelta(exp_time, format='sec')
        t_end = t_begin + dt

        cols = []

        col_e = pyfits.Column(name='ENERGY', format='E', unit='eV',
                              array=self["eobs"].in_units("eV").d)
        col_x = pyfits.Column(name='X', format='D', unit='pixel',
                              array=self["xpix"])
        col_y = pyfits.Column(name='Y', format='D', unit='pixel',
                              array=self["ypix"])

        cols = [col_e, col_x, col_y]

        if "ChannelType" in self.parameters:
             chantype = self.parameters["ChannelType"]
             if chantype == "PHA":
                  cunit = "adu"
             elif chantype == "PI":
                  cunit = "Chan"
             col_ch = pyfits.Column(name=chantype.upper(), format='1J',
                                    unit=cunit, array=self.events[chantype])
             cols.append(col_ch)

             mylog.info("Generating times for events assuming uniform time "
                        "distribution. In future versions this will be made "
                        "more general.")

             time = np.random.uniform(size=self.num_events, low=0.0,
                                      high=float(self.parameters["ExposureTime"]))
             col_t = pyfits.Column(name="TIME", format='1D', unit='s',
                                   array=time)
             cols.append(col_t)

        coldefs = pyfits.ColDefs(cols)
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
        tbhdu.update_ext_name("EVENTS")

        tbhdu.header["MTYPE1"] = "sky"
        tbhdu.header["MFORM1"] = "x,y"
        tbhdu.header["MTYPE2"] = "EQPOS"
        tbhdu.header["MFORM2"] = "RA,DEC"
        tbhdu.header["TCTYP2"] = "RA---TAN"
        tbhdu.header["TCTYP3"] = "DEC--TAN"
        tbhdu.header["TCRVL2"] = float(self.parameters["sky_center"][0])
        tbhdu.header["TCRVL3"] = float(self.parameters["sky_center"][1])
        tbhdu.header["TCDLT2"] = -float(self.parameters["dtheta"])
        tbhdu.header["TCDLT3"] = float(self.parameters["dtheta"])
        tbhdu.header["TCRPX2"] = self.parameters["pix_center"][0]
        tbhdu.header["TCRPX3"] = self.parameters["pix_center"][1]
        tbhdu.header["TLMIN2"] = 0.5
        tbhdu.header["TLMIN3"] = 0.5
        tbhdu.header["TLMAX2"] = 2.*self.parameters["pix_center"][0]-0.5
        tbhdu.header["TLMAX3"] = 2.*self.parameters["pix_center"][1]-0.5
        tbhdu.header["EXPOSURE"] = exp_time
        tbhdu.header["TSTART"] = 0.0
        tbhdu.header["TSTOP"] = exp_time
        tbhdu.header["AREA"] = float(self.parameters["Area"])
        tbhdu.header["D_A"] = float(self.parameters["AngularDiameterDistance"])
        tbhdu.header["REDSHIFT"] = self.parameters["Redshift"]
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header["RADECSYS"] = "FK5"
        tbhdu.header["EQUINOX"] = 2000.0
        tbhdu.header["HDUCLASS"] = "OGIP"
        tbhdu.header["HDUCLAS1"] = "EVENTS"
        tbhdu.header["HDUCLAS2"] = "ACCEPTED"
        tbhdu.header["DATE"] = t_begin.tt.isot
        tbhdu.header["DATE-OBS"] = t_begin.tt.isot
        tbhdu.header["DATE-END"] = t_end.tt.isot
        if "RMF" in self.parameters:
            tbhdu.header["RESPFILE"] = self.parameters["RMF"]
            f = pyfits.open(self.parameters["RMF"])
            nchan = int(f["EBOUNDS"].header["DETCHANS"])
            tbhdu.header["PHA_BINS"] = nchan
            f.close()
        if "ARF" in self.parameters:
            tbhdu.header["ANCRFILE"] = self.parameters["ARF"]
        if "ChannelType" in self.parameters:
            tbhdu.header["CHANTYPE"] = self.parameters["ChannelType"]
        if "Mission" in self.parameters:
            tbhdu.header["MISSION"] = self.parameters["Mission"]
        if "Telescope" in self.parameters:
            tbhdu.header["TELESCOP"] = self.parameters["Telescope"]
        if "Instrument" in self.parameters:
            tbhdu.header["INSTRUME"] = self.parameters["Instrument"]

        hdulist = [pyfits.PrimaryHDU(), tbhdu]

        if "ChannelType" in self.parameters:
            start = pyfits.Column(name='START', format='1D', unit='s',
                                  array=np.array([0.0]))
            stop = pyfits.Column(name='STOP', format='1D', unit='s',
                                 array=np.array([exp_time]))

            tbhdu_gti = pyfits.BinTableHDU.from_columns([start,stop])
            tbhdu_gti.update_ext_name("STDGTI")
            tbhdu_gti.header["TSTART"] = 0.0
            tbhdu_gti.header["TSTOP"] = exp_time
            tbhdu_gti.header["HDUCLASS"] = "OGIP"
            tbhdu_gti.header["HDUCLAS1"] = "GTI"
            tbhdu_gti.header["HDUCLAS2"] = "STANDARD"
            tbhdu_gti.header["RADECSYS"] = "FK5"
            tbhdu_gti.header["EQUINOX"] = 2000.0
            tbhdu_gti.header["DATE"] = t_begin.tt.isot
            tbhdu_gti.header["DATE-OBS"] = t_begin.tt.isot
            tbhdu_gti.header["DATE-END"] = t_end.tt.isot

            hdulist.append(tbhdu_gti)

        pyfits.HDUList(hdulist).writeto(fitsfile, clobber=clobber)

    @parallel_root_only
    def write_simput_file(self, prefix, clobber=False, emin=None, emax=None):
        r"""
        Write events to a SIMPUT file that may be read by the SIMX instrument
        simulator.

        Parameters
        ----------
        prefix : string
            The filename prefix.
        clobber : boolean, optional
            Set to True to overwrite previous files.
        e_min : float, optional
            The minimum energy of the photons to save in keV.
        e_max : float, optional
            The maximum energy of the photons to save in keV.
        """
        pyfits = _astropy.pyfits
        if isinstance(self.parameters["Area"], string_types):
             mylog.error("Writing SIMPUT files is only supported if you didn't convolve with responses.")
             raise TypeError("Writing SIMPUT files is only supported if you didn't convolve with responses.")

        if emin is None:
            emin = self["eobs"].min().value
        if emax is None:
            emax = self["eobs"].max().value

        idxs = np.logical_and(self["eobs"].d >= emin, self["eobs"].d <= emax)
        flux = np.sum(self["eobs"][idxs].in_units("erg")) / \
               self.parameters["ExposureTime"]/self.parameters["Area"]

        col1 = pyfits.Column(name='ENERGY', format='E', array=self["eobs"][idxs].d)
        col2 = pyfits.Column(name='RA', format='D', array=self["xsky"][idxs].d)
        col3 = pyfits.Column(name='DEC', format='D', array=self["ysky"][idxs].d)

        coldefs = pyfits.ColDefs([col1, col2, col3])

        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
        tbhdu.update_ext_name("PHLIST")

        tbhdu.header["HDUCLASS"] = "HEASARC/SIMPUT"
        tbhdu.header["HDUCLAS1"] = "PHOTONS"
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header["EXTVER"] = 1
        tbhdu.header["REFRA"] = 0.0
        tbhdu.header["REFDEC"] = 0.0
        tbhdu.header["TUNIT1"] = "keV"
        tbhdu.header["TUNIT2"] = "deg"
        tbhdu.header["TUNIT3"] = "deg"

        phfile = prefix+"_phlist.fits"

        tbhdu.writeto(phfile, clobber=clobber)

        col1 = pyfits.Column(name='SRC_ID', format='J', array=np.array([1]).astype("int32"))
        col2 = pyfits.Column(name='RA', format='D', array=np.array([0.0]))
        col3 = pyfits.Column(name='DEC', format='D', array=np.array([0.0]))
        col4 = pyfits.Column(name='E_MIN', format='D', array=np.array([float(emin)]))
        col5 = pyfits.Column(name='E_MAX', format='D', array=np.array([float(emax)]))
        col6 = pyfits.Column(name='FLUX', format='D', array=np.array([flux.value]))
        col7 = pyfits.Column(name='SPECTRUM', format='80A', array=np.array([phfile+"[PHLIST,1]"]))
        col8 = pyfits.Column(name='IMAGE', format='80A', array=np.array([phfile+"[PHLIST,1]"]))
        col9 = pyfits.Column(name='SRC_NAME', format='80A', array=np.array(["yt_src"]))

        coldefs = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])

        wrhdu = pyfits.BinTableHDU.from_columns(coldefs)
        wrhdu.update_ext_name("SRC_CAT")

        wrhdu.header["HDUCLASS"] = "HEASARC"
        wrhdu.header["HDUCLAS1"] = "SIMPUT"
        wrhdu.header["HDUCLAS2"] = "SRC_CAT"
        wrhdu.header["HDUVERS"] = "1.1.0"
        wrhdu.header["RADECSYS"] = "FK5"
        wrhdu.header["EQUINOX"] = 2000.0
        wrhdu.header["TUNIT2"] = "deg"
        wrhdu.header["TUNIT3"] = "deg"
        wrhdu.header["TUNIT4"] = "keV"
        wrhdu.header["TUNIT5"] = "keV"
        wrhdu.header["TUNIT6"] = "erg/s/cm**2"

        simputfile = prefix+"_simput.fits"

        wrhdu.writeto(simputfile, clobber=clobber)

    @parallel_root_only
    def write_h5_file(self, h5file):
        """
        Write an :class:`~pyxsim.event_list.EventList` to the HDF5 file given by *h5file*.
        """
        f = h5py.File(h5file, "w")

        p = f.create_group("parameters")
        p.create_dataset("exp_time", data=float(self.parameters["ExposureTime"]))
        p.create_dataset("area", data=float(self.parameters["Area"]))
        p.create_dataset("redshift", data=self.parameters["Redshift"])
        p.create_dataset("d_a", data=float(self.parameters["AngularDiameterDistance"]))
        if "ARF" in self.parameters:
            p.create_dataset("arf", data=self.parameters["ARF"])
        if "RMF" in self.parameters:
            p.create_dataset("rmf", data=self.parameters["RMF"])
        if "ChannelType" in self.parameters:
            p.create_dataset("channel_type", data=self.parameters["ChannelType"])
        if "Mission" in self.parameters:
            p.create_dataset("mission", data=self.parameters["Mission"])
        if "Telescope" in self.parameters:
            p.create_dataset("telescope", data=self.parameters["Telescope"])
        if "Instrument" in self.parameters:
            p.create_dataset("instrument", data=self.parameters["Instrument"])
        p.create_dataset("sky_center", data=self.parameters["sky_center"].d)
        p.create_dataset("pix_center", data=self.parameters["pix_center"])
        p.create_dataset("dtheta", data=float(self.parameters["dtheta"]))

        d = f.create_group("data")
        d.create_dataset("xpix", data=self["xpix"])
        d.create_dataset("ypix", data=self["ypix"])
        d.create_dataset("xsky", data=self["xsky"].d)
        d.create_dataset("ysky", data=self["ysky"].d)
        d.create_dataset("eobs", data=self["eobs"].d)
        if "PI" in self.events:
            d.create_dataset("pi", data=self.events["PI"])
        if "PHA" in self.events:
            d.create_dataset("pha", data=self.events["PHA"])

        f.close()

    @parallel_root_only
    def write_fits_image(self, imagefile, clobber=False,
                         emin=None, emax=None):
        r"""
        Generate a image by binning X-ray counts and write it to a FITS file.

        Parameters
        ----------
        imagefile : string
            The name of the image file to write.
        clobber : boolean, optional
            Set to True to overwrite a previous file.
        emin : float, optional
            The minimum energy of the photons to put in the image, in keV.
        emax : float, optional
            The maximum energy of the photons to put in the image, in keV.
        """
        if emin is None:
            mask_emin = np.ones(self.num_events, dtype='bool')
        else:
            mask_emin = self["eobs"].d > emin
        if emax is None:
            mask_emax = np.ones(self.num_events, dtype='bool')
        else:
            mask_emax = self["eobs"].d < emax

        mask = np.logical_and(mask_emin, mask_emax)

        nx = int(2*self.parameters["pix_center"][0]-1.)
        ny = int(2*self.parameters["pix_center"][1]-1.)

        xbins = np.linspace(0.5, float(nx)+0.5, nx+1, endpoint=True)
        ybins = np.linspace(0.5, float(ny)+0.5, ny+1, endpoint=True)

        H, xedges, yedges = np.histogram2d(self["xpix"][mask],
                                           self["ypix"][mask],
                                           bins=[xbins,ybins])

        hdu = _astropy.pyfits.PrimaryHDU(H.T)

        hdu.header["MTYPE1"] = "EQPOS"
        hdu.header["MFORM1"] = "RA,DEC"
        hdu.header["CTYPE1"] = "RA---TAN"
        hdu.header["CTYPE2"] = "DEC--TAN"
        hdu.header["CRPIX1"] = 0.5*(nx+1)
        hdu.header["CRPIX2"] = 0.5*(nx+1)
        hdu.header["CRVAL1"] = float(self.parameters["sky_center"][0])
        hdu.header["CRVAL2"] = float(self.parameters["sky_center"][1])
        hdu.header["CUNIT1"] = "deg"
        hdu.header["CUNIT2"] = "deg"
        hdu.header["CDELT1"] = -float(self.parameters["dtheta"])
        hdu.header["CDELT2"] = float(self.parameters["dtheta"])
        hdu.header["EXPOSURE"] = float(self.parameters["ExposureTime"])

        hdu.writeto(imagefile, clobber=clobber)

    @parallel_root_only
    def write_spectrum(self, specfile, bin_type="channel", emin=0.1,
                       emax=10.0, nchan=2000, clobber=False):
        r"""
        Bin event energies into a spectrum and write it to a FITS binary table. Can bin
        on energy or channel. In that case, the spectral binning will be determined by 
        the RMF binning.

        Parameters
        ----------
        specfile : string
            The name of the FITS file to be written.
        bin_type : string, optional
            Bin on "energy" or "channel". If an RMF is detected, channel information will be
            imported from it. 
        emin : float, optional
            The minimum energy of the spectral bins in keV. Only used if binning without an RMF.
        emax : float, optional
            The maximum energy of the spectral bins in keV. Only used if binning without an RMF.
        nchan : integer, optional
            The number of channels. Only used if binning without an RMF.
        """
        pyfits = _astropy.pyfits
        if bin_type == "channel" and "ChannelType" in self.parameters:
            spectype = self.parameters["ChannelType"]
            f = pyfits.open(self.parameters["RMF"])
            nchan = int(f["EBOUNDS"].header["DETCHANS"])
            num = 0
            for i in range(1,len(f["EBOUNDS"].columns)+1):
                if f["EBOUNDS"].header["TTYPE%d" % i] == "CHANNEL":
                    num = i
                    break
            if num > 0:
                tlmin = "TLMIN%d" % num
                cmin = int(f["EBOUNDS"].header[tlmin])
            else:
                mylog.warning("Cannot determine minimum allowed value for channel. " +
                              "Setting to 0, which may be wrong.")
                cmin = 0
            f.close()
            minlength = nchan
            if cmin == 1: minlength += 1
            spec = np.bincount(self[spectype],minlength=minlength)
            if cmin == 1: spec = spec[1:]
            bins = (np.arange(nchan)+cmin).astype("int32")
        else:
            espec = self["eobs"].d
            erange = (emin, emax)
            spec, ee = np.histogram(espec, bins=nchan, range=erange)
            if bin_type == "energy":
                bins = 0.5*(ee[1:]+ee[:-1])
                spectype = "energy"
            else:
                mylog.info("Events haven't been convolved with an RMF, so assuming "
                           "a perfect response and %d PI channels." % nchan)
                bins = (np.arange(nchan)+1).astype("int32")
                spectype = "pi"

        col1 = pyfits.Column(name='CHANNEL', format='1J', array=bins)
        col2 = pyfits.Column(name=spectype.upper(), format='1D', array=bins.astype("float64"))
        col3 = pyfits.Column(name='COUNTS', format='1J', array=spec.astype("int32"))
        col4 = pyfits.Column(name='COUNT_RATE', format='1D', array=spec/float(self.parameters["ExposureTime"]))

        coldefs = pyfits.ColDefs([col1, col2, col3, col4])

        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
        tbhdu.update_ext_name("SPECTRUM")

        tbhdu.header["DETCHANS"] = spec.shape[0]
        tbhdu.header["TOTCTS"] = spec.sum()
        tbhdu.header["EXPOSURE"] = float(self.parameters["ExposureTime"])
        tbhdu.header["LIVETIME"] = float(self.parameters["ExposureTime"])
        tbhdu.header["CONTENT"] = spectype
        tbhdu.header["HDUCLASS"] = "OGIP"
        tbhdu.header["HDUCLAS1"] = "SPECTRUM"
        tbhdu.header["HDUCLAS2"] = "TOTAL"
        tbhdu.header["HDUCLAS3"] = "TYPE:I"
        tbhdu.header["HDUCLAS4"] = "COUNT"
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header["HDUVERS1"] = "1.1.0"
        tbhdu.header["CHANTYPE"] = spectype
        tbhdu.header["BACKFILE"] = "none"
        tbhdu.header["CORRFILE"] = "none"
        tbhdu.header["POISSERR"] = True
        if "RMF" in self.parameters:
            tbhdu.header["RESPFILE"] = self.parameters["RMF"]
        else:
            tbhdu.header["RESPFILE"] = "none"
        if "ARF" in self.parameters:
            tbhdu.header["ANCRFILE"] = self.parameters["ARF"]
        else:
            tbhdu.header["ANCRFILE"] = "none"
        if "Mission" in self.parameters:
            tbhdu.header["MISSION"] = self.parameters["Mission"]
        else:
            tbhdu.header["MISSION"] = "none"
        if "Telescope" in self.parameters:
            tbhdu.header["TELESCOP"] = self.parameters["Telescope"]
        else:
            tbhdu.header["TELESCOP"] = "none"
        if "Instrument" in self.parameters:
            tbhdu.header["INSTRUME"] = self.parameters["Instrument"]
        else:
            tbhdu.header["INSTRUME"] = "none"
        tbhdu.header["AREASCAL"] = 1.0
        tbhdu.header["CORRSCAL"] = 0.0
        tbhdu.header["BACKSCAL"] = 1.0

        hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])

        hdulist.writeto(specfile, clobber=clobber)