"""
Classes for generating lists of detected events
"""
import numpy as np
from pyxsim.utils import mylog
from astropy.io import fits
import astropy.wcs as pywcs
import h5py
from pyxsim.utils import parse_value
from soxs.simput import write_photon_list


class EventList(object):

    def __init__(self, filespec):
        from glob import glob
        if filespec.endswith(".h5"):
            if "*" in filespec:
                filenames = glob(filespec)
            else:
                filenames = [filespec]
        elif isinstance(filespec, list):
            if not np.all([fn.endswith("*.h5") for fn in filespec]):
                raise RuntimeError("Not all filenames are valid!")
            filenames = filespec
        else:
            raise RuntimeError(f"Invalid EventList file spec: {filespec}")
        self.filenames = filenames
        self.parameters = {}
        self.num_events = []
        for i, fn in enumerate(self.filenames):
            with h5py.File(fn, "r") as f:
                p = f["parameters"]
                self.num_events.append(f["data"]["xsky"].size)
                if i == 0:
                    for field in p:
                        self.parameters[field] = p[field][()]
        self.tot_num_events = np.sum(self.num_events)

    def write_fits_file(self, fitsfile, fov, nx, overwrite=False):
        """
        Write events to a FITS binary table file. The result is a
        standard "event file" which can be processed by standard
        X-ray analysis tools.

        Parameters
        ----------
        fitsfile : string
            The name of the event file to write.
        fov : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The field of view of the event file. If units are not 
            provided, they are assumed to be in arcminutes.
        nx : integer
            The resolution of the image (number of pixels on a side). 
        overwrite : boolean, optional
            Set to True to overwrite a previous file.
        """
        from astropy.time import Time, TimeDelta

        fov = parse_value(fov, "arcmin")

        t_begin = Time.now()
        dt = TimeDelta(self.parameters["exp_time"], format='sec')
        t_end = t_begin + dt

        dtheta = fov.to("deg").v / nx

        wcs = pywcs.WCS(naxis=2)
        wcs.wcs.crpix = [0.5*(nx+1)]*2
        wcs.wcs.crval = self.parameters["sky_center"]
        wcs.wcs.cdelt = [-dtheta, dtheta]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.cunit = ["deg"]*2

        xx, yy = wcs.wcs_world2pix(self["xsky"].d, self["ysky"].d, 1)

        keepx = np.logical_and(xx >= 0.5, xx <= float(nx)+0.5)
        keepy = np.logical_and(yy >= 0.5, yy <= float(nx)+0.5)
        keep = np.logical_and(keepx, keepy)

        n_events = keep.sum()

        mylog.info("Threw out %d events because " % (xx.size-n_events) +
                   "they fell outside the field of view.")

        col_e = fits.Column(name='ENERGY', format='E', unit='eV',
                            array=events["eobs"].in_units("eV").d[keep])
        col_x = fits.Column(name='X', format='D', unit='pixel',
                            array=xx[keep])
        col_y = fits.Column(name='Y', format='D', unit='pixel',
                            array=yy[keep])

        cols = [col_e, col_x, col_y]

        coldefs = fits.ColDefs(cols)
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = "EVENTS"

        tbhdu.header["MTYPE1"] = "sky"
        tbhdu.header["MFORM1"] = "x,y"
        tbhdu.header["MTYPE2"] = "EQPOS"
        tbhdu.header["MFORM2"] = "RA,DEC"
        tbhdu.header["TCTYP2"] = "RA---TAN"
        tbhdu.header["TCTYP3"] = "DEC--TAN"
        tbhdu.header["TCRVL2"] = self.parameters["sky_center"][0]
        tbhdu.header["TCRVL3"] = self.parameters["sky_center"][1]
        tbhdu.header["TCDLT2"] = -dtheta
        tbhdu.header["TCDLT3"] = dtheta
        tbhdu.header["TCRPX2"] = 0.5*(nx+1)
        tbhdu.header["TCRPX3"] = 0.5*(nx+1)
        tbhdu.header["TLMIN2"] = 0.5
        tbhdu.header["TLMIN3"] = 0.5
        tbhdu.header["TLMAX2"] = float(nx)+0.5
        tbhdu.header["TLMAX3"] = float(nx)+0.5
        tbhdu.header["EXPOSURE"] = self.parameters["exp_time"]
        tbhdu.header["TSTART"] = 0.0
        tbhdu.header["TSTOP"] = self.parameters["exp_time"]
        tbhdu.header["AREA"] = self.parameters["area"]
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header["RADECSYS"] = "FK5"
        tbhdu.header["EQUINOX"] = 2000.0
        tbhdu.header["HDUCLASS"] = "OGIP"
        tbhdu.header["HDUCLAS1"] = "EVENTS"
        tbhdu.header["HDUCLAS2"] = "ACCEPTED"
        tbhdu.header["DATE"] = t_begin.tt.isot
        tbhdu.header["DATE-OBS"] = t_begin.tt.isot
        tbhdu.header["DATE-END"] = t_end.tt.isot

        hdulist = [fits.PrimaryHDU(), tbhdu]

        fits.HDUList(hdulist).writeto(fitsfile, overwrite=overwrite)

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

        mylog.info(f"Writing SIMPUT catalog file {prefix}_simput.fits "
                   f"and SIMPUT photon list file {prefix}_phlist.fits.")

        if emin is None and emax is None:
            idxs = slice(None, None, None)
        else:
            if emin is None:
                emin = events["eobs"].min().value
            if emax is None:
                emax = events["eobs"].max().value
            idxs = np.logical_and(events["eobs"].d >= emin, events["eobs"].d <= emax)

        flux = np.sum(events["eobs"][idxs]).to("erg") / \
               self.parameters["exp_time"]/self.parameters["area"]

        write_photon_list(prefix, prefix, flux.v, events["xsky"][idxs].d,
                          events["ysky"][idxs].d, events["eobs"][idxs].d,
                          overwrite=overwrite)

    def write_fits_image(self, imagefile, fov, nx, emin=None, 
                         emax=None, overwrite=False):
        r"""
        Generate a image by binning X-ray counts and write it to a FITS file.

        Parameters
        ----------
        imagefile : string
            The name of the image file to write.
        fov : float, (value, unit) tuple, :class:`~yt.units.yt_array.YTQuantity`, or :class:`~astropy.units.Quantity`
            The field of view of the image. If units are not provided, they
            are assumed to be in arcminutes.
        nx : integer
            The resolution of the image (number of pixels on a side). 
        emin : float, optional
            The minimum energy of the photons to put in the image, in keV.
        emax : float, optional
            The maximum energy of the photons to put in the image, in keV.
        overwrite : boolean, optional
            Set to True to overwrite a previous file.
        """
        fov = parse_value(fov, "arcmin")

        if emin is None:
            mask_emin = np.ones(self.num_events, dtype='bool')
        else:
            mask_emin = self["eobs"].d > emin
        if emax is None:
            mask_emax = np.ones(self.num_events, dtype='bool')
        else:
            mask_emax = self["eobs"].d < emax

        mask = np.logical_and(mask_emin, mask_emax)

        dtheta = fov.to_value("deg").v/nx

        xbins = np.linspace(0.5, float(nx)+0.5, nx+1, endpoint=True)
        ybins = np.linspace(0.5, float(nx)+0.5, nx+1, endpoint=True)

        wcs = pywcs.WCS(naxis=2)
        wcs.wcs.crpix = [0.5*(nx+1)]*2
        wcs.wcs.crval = self.parameters["sky_center"].d
        wcs.wcs.cdelt = [-dtheta, dtheta]
        wcs.wcs.ctype = ["RA---TAN","DEC--TAN"]
        wcs.wcs.cunit = ["deg"]*2

        xx, yy = wcs.wcs_world2pix(self["xsky"].d, self["ysky"].d, 1)

        H, xedges, yedges = np.histogram2d(xx[mask], yy[mask],
                                           bins=[xbins, ybins])

        hdu = fits.PrimaryHDU(H.T)

        hdu.header["MTYPE1"] = "EQPOS"
        hdu.header["MFORM1"] = "RA,DEC"
        hdu.header["CTYPE1"] = "RA---TAN"
        hdu.header["CTYPE2"] = "DEC--TAN"
        hdu.header["CRPIX1"] = 0.5*(nx+1)
        hdu.header["CRPIX2"] = 0.5*(nx+1)
        hdu.header["CRVAL1"] = self.parameters["sky_center"][0]
        hdu.header["CRVAL2"] = self.parameters["sky_center"][1]
        hdu.header["CUNIT1"] = "deg"
        hdu.header["CUNIT2"] = "deg"
        hdu.header["CDELT1"] = -dtheta
        hdu.header["CDELT2"] = dtheta
        hdu.header["EXPOSURE"] = self.parameters["exp_time"]

        hdu.writeto(imagefile, overwrite=overwrite)

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
        espec = self["eobs"]
        spec, ee = np.histogram(espec, bins=nchan, range=(emin, emax))
        bins = 0.5*(ee[1:]+ee[:-1])

        col1 = fits.Column(name='CHANNEL', format='1J', 
                           array=np.arange(nchan).astype('int32')+1)
        col2 = fits.Column(name='ENERGY', format='1D', 
                           array=bins.astype("float64"))
        col3 = fits.Column(name='COUNTS', format='1J', 
                           array=spec.astype("int32"))
        col4 = fits.Column(name='COUNT_RATE', format='1D', 
                           array=spec/self.parameters["exp_time"])

        coldefs = fits.ColDefs([col1, col2, col3, col4])

        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = "SPECTRUM"

        tbhdu.header["DETCHANS"] = spec.shape[0]
        tbhdu.header["TOTCTS"] = spec.sum()
        tbhdu.header["EXPOSURE"] = self.parameters["exp_time"]
        tbhdu.header["LIVETIME"] = self.parameters["exp_time"]
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
