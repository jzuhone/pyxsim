"""
Classes for generating lists of detected events
"""
from six import string_types
import numpy as np
from yt.funcs import ensure_list
from pyxsim.utils import mylog
try:
    from yt.utilities.fits_image import assert_same_wcs
except ImportError:
    from yt.visualization.fits_image import assert_same_wcs
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.units.yt_array import YTQuantity, YTArray, uconcatenate
from yt.utilities.on_demand_imports import _astropy
import h5py
from pyxsim.utils import force_unicode, validate_parameters, parse_value
from pyxsim.responses import RedistributionMatrixFile
import os

class EventList(object):

    def __init__(self, events, parameters, wcs=None):
        self.events = events
        self.parameters = parameters
        self.num_events = events["xpix"].shape[0]
        if wcs is None:
            self.wcs = _astropy.pywcs.WCS(naxis=2)
            self.wcs.wcs.crpix = parameters["pix_center"]
            self.wcs.wcs.crval = parameters["sky_center"].d
            self.wcs.wcs.cdelt = [-parameters["dtheta"].value, parameters["dtheta"].value]
            self.wcs.wcs.ctype = ["RA---TAN","DEC--TAN"]
            self.wcs.wcs.cunit = ["deg"]*2
        else:
            self.wcs = wcs

    @classmethod
    def create_empty_list(cls, exp_time, area, wcs, parameters=None):
        events = {"xpix": np.array([]),
                  "ypix": np.array([]),
                  "eobs": YTArray([], "keV")}
        if parameters is None:
            parameters = {}
        parameters["ExposureTime"] = parse_value(exp_time, "s")
        parameters["Area"] = parse_value(area, "cm**2")
        parameters["pix_center"] = wcs.wcs.crpix[:]
        parameters["sky_center"] = YTArray(wcs.wcs.crval[:], "deg")
        parameters["dtheta"] = YTQuantity(wcs.wcs.cdelt[0], "deg")
        return cls(events, parameters)

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

    def add_point_sources(self, positions, energy_bins, spectra,
                          prng=None, absorb_model=None):
        r"""
        Add point source events to an :class:`~pyxsim.event_list.EventList`.
        Returns a new :class:`~pyxsim.event_list.EventList`.

        Parameters
        ----------
        positions : array of source positions, shape 2xN
            The positions of the point sources in RA, Dec, where N is the
            number of point sources. Coordinates should be in degrees.
        energy_bins : :class:`~yt.units.yt_array.YTArray` with units of keV, shape M+1
            The edges of the energy bins for the spectra, where M is the number of
            bins
        spectra : list (size N) of :class:`~yt.units.yt_array.YTArray`\s with units of photons/s/cm**2, each with shape M
            The spectra for the point sources, where M is the number of bins and N is
            the number of point sources
        prng : :class:`~numpy.random.RandomState` object or :mod:`numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        absorb_model : :class:`~pyxsim.spectral_models.AbsorptionModel` 
            A model for foreground galactic absorption.
        """
        if prng is None:
            prng = np.random

        spectra = ensure_list(spectra)
        positions = ensure_list(positions)

        x = [self.events["xpix"]]
        y = [self.events["ypix"]]
        e = [self.events["eobs"]]

        for pos, spectrum in zip(positions, spectra):
            eobs = self._add_events(energy_bins, spectrum, prng, absorb_model)
            xpix, ypix = self.wcs.wcs_world2pix(pos[0], pos[1], 1)
            ne = len(eobs)
            x.append([xpix]*ne)
            y.append([ypix]*ne)
            e.append(eobs)

        events = {}
        events["xpix"] = uconcatenate(x)
        events["ypix"] = uconcatenate(y)
        events["eobs"] = uconcatenate(e)

        return EventList(events, self.parameters)

    def add_background(self, energy_bins, spectrum,
                       prng=None, absorb_model=None):
        r"""
        Add background events to an :class:`~pyxsim.event_list.EventList`.
        Returns a new :class:`~pyxsim.event_list.EventList`.

        Parameters
        ----------
        energy_bins : :class:`~yt.units.yt_array.YTArray` with units of keV, size M+1
            The edges of the energy bins for the spectra, where M is the number of
            bins
        spectrum : :class:`~yt.units.yt_array.YTArray` with units of photons/s/cm**2, size M
            The spectrum for the background, where M is the number of bins.
        prng : :class:`~numpy.random.RandomState` object or :mod:`numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        absorb_model : :class:`~pyxsim.spectral_models.AbsorptionModel` 
            A model for foreground galactic absorption.
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
        exp_time = self.parameters["ExposureTime"]
        area = self.parameters["Area"]
        flux = spectrum.sum()
        num_photons = np.uint64(exp_time*area*flux)
        cumspec = np.cumsum(spectrum)
        cumspec = np.insert(cumspec, 0, 0.0)
        cumspec /= cumspec[-1]
        randvec = prng.uniform(size=num_photons)
        randvec.sort()
        e = YTArray(np.interp(randvec, cumspec, ebins), "keV")

        if absorb_model is None:
            detected = np.ones(e.shape, dtype='bool')
        else:
            detected = absorb_model.absorb_photons(e, prng=prng)

        mylog.info("Adding %d new events." % detected.sum())

        return e[detected]

    def filter_events(self, region):
        """
        Filter events using a ds9 *region*. Requires the `pyregion <http://pyregion.readthedocs.org/en/latest/>`_ package.
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
        if "redshift" in p:
            parameters["Redshift"] = p["redshift"].value
        if "d_a" in p:
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
        if "REDSHIFT" in tblhdu.header:
            parameters["Redshift"] = tblhdu.header["REDSHIFT"]
        if "D_A" in tblhdu.header:
            parameters["AngularDiameterDistance"] = YTQuantity(tblhdu.header["D_A"], "Mpc")
        if "RMF" in tblhdu.header:
            parameters["RMF"] = tblhdu.header["RMF"]
        if "ARF" in tblhdu.header:
            parameters["ARF"] = tblhdu.header["ARF"]
        if "CHANTYPE" in tblhdu.header:
            parameters["ChannelType"] = tblhdu.header["CHANTYPE"]
        if "MISSION" in tblhdu.header:
            parameters["Mission"] = tblhdu.header["MISSION"]
        if "TELESCOP" in tblhdu.header:
            parameters["Telescope"] = tblhdu.header["TELESCOP"]
        if "INSTRUME" in tblhdu.header:
            parameters["Instrument"] = tblhdu.header["INSTRUME"]
        parameters["sky_center"] = YTArray([tblhdu.header["TCRVL2"], tblhdu.header["TCRVL3"]], "deg")
        parameters["pix_center"] = np.array([tblhdu.header["TCRVL2"], tblhdu.header["TCRVL3"]])
        parameters["dtheta"] = YTQuantity(tblhdu.header["TCRVL3"], "deg")
        events["xpix"] = tblhdu.data["X"]
        events["ypix"] = tblhdu.data["Y"]
        events["eobs"] = YTArray(tblhdu.data["ENERGY"]/1000., "keV")
        if "PI" in tblhdu.columns.names:
            events["PI"] = tblhdu.data["PI"]
        if "PHA" in tblhdu.columns.names:
            events["PHA"] = tblhdu.data["PHA"]

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
        tbhdu.name = "EVENTS"

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
        if "ChannelType" in self.parameters:
            rmf = RedistributionMatrixFile(self.parameters["RMF"])
            tbhdu.header["TLMIN4"] = rmf.cmin
            tbhdu.header["TLMAX4"] = rmf.cmax
            tbhdu.header["RESPFILE"] = os.path.split(self.parameters["RMF"])[-1]
            tbhdu.header["PHA_BINS"] = rmf.n_ch
        tbhdu.header["EXPOSURE"] = exp_time
        tbhdu.header["TSTART"] = 0.0
        tbhdu.header["TSTOP"] = exp_time
        tbhdu.header["AREA"] = float(self.parameters["Area"])
        if "AngularDiameterDistance" in self.parameters:
            tbhdu.header["D_A"] = float(self.parameters["AngularDiameterDistance"])
        if "Redshift" in self.parameters:
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
        if "ARF" in self.parameters:
            tbhdu.header["ANCRFILE"] = os.path.split(self.parameters["ARF"])[-1]
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
            tbhdu_gti.name = "STDGTI"
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
        tbhdu.name = "PHLIST"

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
        wrhdu.name = "SRC_CAT"

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
        if "Redshift" in self.parameters:
            p.create_dataset("redshift", data=self.parameters["Redshift"])
        if "AngularDiameterDistance" in self.parameters:
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
        on energy or channel. In the latter case, the spectral binning will be determined by
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
            rmf = RedistributionMatrixFile(self.parameters["RMF"])
            minlength = rmf.n_ch
            if rmf.cmin == 1: minlength += 1
            spec = np.bincount(self[spectype], minlength=minlength)
            if rmf.cmin == 1: spec = spec[1:]
            bins = (np.arange(rmf.n_ch)+rmf.cmin).astype("int32")
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
        tbhdu.name = "SPECTRUM"

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
            tbhdu.header["RESPFILE"] = os.path.split(self.parameters["RMF"])[-1]
        else:
            tbhdu.header["RESPFILE"] = "none"
        if "ARF" in self.parameters:
            tbhdu.header["ANCRFILE"] = os.path.split(self.parameters["ARF"])[-1]
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
