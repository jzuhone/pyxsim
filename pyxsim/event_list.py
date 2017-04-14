"""
Classes for generating lists of detected events
"""
import numpy as np
from yt.funcs import ensure_list
from pyxsim.utils import mylog
try:
    from yt.visualization.fits_image import assert_same_wcs
except ImportError:
    from yt.utilities.fits_image import assert_same_wcs
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.units.yt_array import YTQuantity, YTArray, uconcatenate
from yt.utilities.on_demand_imports import _astropy
import h5py
from pyxsim.utils import validate_parameters, parse_value, \
    ParameterDict
from soxs.simput import write_photon_list

old_parameter_keys = {"ExposureTime": "exp_time",
                      "Area": "area",
                      "Redshift": "redshift",
                      "AngularDiameterDistance": "d_a"}

class EventList(object):

    def __init__(self, events, parameters, wcs=None):
        self.events = events
        self.parameters = ParameterDict(parameters, "EventList", old_parameter_keys)
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
        parameters["exp_time"] = parse_value(exp_time, "s")
        parameters["area"] = parse_value(area, "cm**2")
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
        parameters["exp_time"] = YTQuantity(p["exp_time"].value, "s")
        parameters["area"] = YTQuantity(p["area"].value, "cm**2")
        if "redshift" in p:
            parameters["redshift"] = p["redshift"].value
        if "d_a" in p:
            parameters["d_a"] = YTQuantity(p["d_a"].value, "Mpc")
        parameters["sky_center"] = YTArray(p["sky_center"][:], "deg")
        parameters["dtheta"] = YTQuantity(p["dtheta"].value, "deg")
        parameters["pix_center"] = p["pix_center"][:]

        d = f["/data"]
        events["xpix"] = d["xpix"][:]
        events["ypix"] = d["ypix"][:]
        events["eobs"] = YTArray(d["eobs"][:], "keV")

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

        parameters["exp_time"] = YTQuantity(tblhdu.header["EXPOSURE"], "s")
        parameters["area"] = YTQuantity(tblhdu.header["AREA"], "cm**2")
        if "REDSHIFT" in tblhdu.header:
            parameters["redshift"] = tblhdu.header["REDSHIFT"]
        if "D_A" in tblhdu.header:
            parameters["d_a"] = YTQuantity(tblhdu.header["D_A"], "Mpc")
        parameters["sky_center"] = YTArray([tblhdu.header["TCRVL2"], tblhdu.header["TCRVL3"]], "deg")
        parameters["pix_center"] = np.array([tblhdu.header["TCRVL2"], tblhdu.header["TCRVL3"]])
        parameters["dtheta"] = YTQuantity(tblhdu.header["TCRVL3"], "deg")
        events["xpix"] = tblhdu.data["X"]
        events["ypix"] = tblhdu.data["Y"]
        events["eobs"] = YTArray(tblhdu.data["ENERGY"]/1000., "keV")

        return cls(events, parameters)

    @parallel_root_only
    def write_fits_file(self, fitsfile, overwrite=False):
        """
        Write events to a FITS binary table file with filename *fitsfile*.
        Set *overwrite* to True if you need to overwrite a previous file.
        """
        pyfits = _astropy.pyfits

        exp_time = float(self.parameters["exp_time"])

        col_e = pyfits.Column(name='ENERGY', format='E', unit='eV',
                              array=self["eobs"].in_units("eV").d)
        col_x = pyfits.Column(name='X', format='D', unit='pixel',
                              array=self["xpix"])
        col_y = pyfits.Column(name='Y', format='D', unit='pixel',
                              array=self["ypix"])

        cols = [col_e, col_x, col_y]

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
        tbhdu.header["EXPOSURE"] = exp_time
        tbhdu.header["TSTART"] = 0.0
        tbhdu.header["TSTOP"] = exp_time
        tbhdu.header["AREA"] = float(self.parameters["area"])
        if "d_a" in self.parameters:
            tbhdu.header["D_A"] = float(self.parameters["d_a"])
        if "redshift" in self.parameters:
            tbhdu.header["REDSHIFT"] = self.parameters["redshift"]
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header["RADECSYS"] = "FK5"
        tbhdu.header["EQUINOX"] = 2000.0
        tbhdu.header["HDUCLASS"] = "OGIP"
        tbhdu.header["HDUCLAS1"] = "EVENTS"
        tbhdu.header["HDUCLAS2"] = "ACCEPTED"

        hdulist = [pyfits.PrimaryHDU(), tbhdu]

        pyfits.HDUList(hdulist).writeto(fitsfile, overwrite=overwrite)

    @parallel_root_only
    def write_simput_file(self, prefix, overwrite=False, emin=None, emax=None):
        r"""
        Write events to a SIMPUT file that may be read by the SIMX instrument
        simulator.

        Parameters
        ----------
        prefix : string
            The filename prefix.
        overwrite : boolean, optional
            Set to True to overwrite previous files.
        e_min : float, optional
            The minimum energy of the photons to save in keV.
        e_max : float, optional
            The maximum energy of the photons to save in keV.
        """
        if emin is None:
            emin = self["eobs"].min().value
        if emax is None:
            emax = self["eobs"].max().value

        idxs = np.logical_and(self["eobs"].d >= emin, self["eobs"].d <= emax)
        flux = np.sum(self["eobs"][idxs].in_units("erg")) / \
               self.parameters["exp_time"]/self.parameters["area"]

        write_photon_list(prefix, prefix, flux.v, self["xsky"][idxs].d, self["ysky"][idxs].d,
                          self["eobs"][idxs].d, overwrite=overwrite)

    @parallel_root_only
    def write_h5_file(self, h5file):
        """
        Write an :class:`~pyxsim.event_list.EventList` to the HDF5 file given by *h5file*.
        """
        f = h5py.File(h5file, "w")

        p = f.create_group("parameters")
        p.create_dataset("exp_time", data=float(self.parameters["exp_time"]))
        p.create_dataset("area", data=float(self.parameters["area"]))
        if "redshift" in self.parameters:
            p.create_dataset("redshift", data=self.parameters["redshift"])
        if "d_a" in self.parameters:
            p.create_dataset("d_a", data=float(self.parameters["d_a"]))
        p.create_dataset("sky_center", data=self.parameters["sky_center"].d)
        p.create_dataset("pix_center", data=self.parameters["pix_center"])
        p.create_dataset("dtheta", data=float(self.parameters["dtheta"]))

        d = f.create_group("data")
        d.create_dataset("xpix", data=self["xpix"])
        d.create_dataset("ypix", data=self["ypix"])
        d.create_dataset("xsky", data=self["xsky"].d)
        d.create_dataset("ysky", data=self["ysky"].d)
        d.create_dataset("eobs", data=self["eobs"].d)

        f.close()

    @parallel_root_only
    def write_fits_image(self, imagefile, overwrite=False,
                         emin=None, emax=None):
        r"""
        Generate a image by binning X-ray counts and write it to a FITS file.

        Parameters
        ----------
        imagefile : string
            The name of the image file to write.
        overwrite : boolean, optional
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
        hdu.header["EXPOSURE"] = float(self.parameters["exp_time"])

        hdu.writeto(imagefile, overwrite=overwrite)

    @parallel_root_only
    def write_spectrum(self, specfile, emin, emax, nchan, overwrite=False):
        r"""
        Bin event energies into a spectrum and write it to a FITS binary table. 
        This is for an *unconvolved* spectrum.

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
        pyfits = _astropy.pyfits
        espec = self["eobs"].d
        spec, ee = np.histogram(espec, bins=nchan, range=(emin, emax))
        bins = 0.5*(ee[1:]+ee[:-1])

        col1 = pyfits.Column(name='CHANNEL', format='1J', array=np.arange(nchan).astype('int32')+1)
        col2 = pyfits.Column(name='ENERGY', format='1D', array=bins.astype("float64"))
        col3 = pyfits.Column(name='COUNTS', format='1J', array=spec.astype("int32"))
        col4 = pyfits.Column(name='COUNT_RATE', format='1D', array=spec/float(self.parameters["exp_time"]))

        coldefs = pyfits.ColDefs([col1, col2, col3, col4])

        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = "SPECTRUM"

        tbhdu.header["DETCHANS"] = spec.shape[0]
        tbhdu.header["TOTCTS"] = spec.sum()
        tbhdu.header["EXPOSURE"] = float(self.parameters["exp_time"])
        tbhdu.header["LIVETIME"] = float(self.parameters["exp_time"])
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

        hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])

        hdulist.writeto(specfile, overwrite=overwrite)
