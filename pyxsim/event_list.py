"""
Classes for generating lists of detected events
"""
import os

import astropy.wcs as pywcs
import h5py
import numpy as np
from astropy.io import fits

from pyxsim.utils import mylog, parse_value


class EventList:
    def __init__(self, filespec):
        """
        Read an EventList from disk so it can be exported to
        other formats.

        Parameters
        ----------
        filespec : str
            A filename, list of filenames, or globbable path
            that yields a single or list of HDF5 files containing
            event data.
        """
        from glob import glob

        if filespec.endswith(".h5"):
            filenames = glob(filespec)
        elif isinstance(filespec, list):
            if not np.all([fn.endswith(".h5") for fn in filespec]):
                raise RuntimeError("Not all filenames are valid!")
            filenames = filespec
        else:
            raise RuntimeError(f"Invalid EventList file spec: {filespec}")
        self.filenames = filenames
        self.filenames.sort()
        self.parameters = {}
        self.info = {}
        self.num_events = []
        for i, fn in enumerate(self.filenames):
            with h5py.File(fn, "r") as f:
                p = f["parameters"]
                info = f["info"]
                self.num_events.append(f["data"]["xsky"].size)
                if i == 0:
                    for field in p:
                        if isinstance(p[field][()], (str, bytes)):
                            self.parameters[field] = p[field].asstr()[()]
                        else:
                            self.parameters[field] = p[field][()]
                    for k, v in info.attrs.items():
                        self.info[k] = v
        self.tot_num_events = np.sum(self.num_events)
        self.observer = self.parameters.get("observer", "external")

    def write_fits_file(self, fitsfile, fov, nx, overwrite=False):
        """
        Write events to a FITS binary table file. The result is an
        "event file" which can be opened in ds9, binned into a
        spectrum, etc., but NOTE that this does NOT represent an
        actual observation since no instrumental effects are
        included.

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
        dt = TimeDelta(self.parameters["exp_time"], format="sec")
        t_end = t_begin + dt

        dtheta = fov.to_value("deg") / nx

        wcs = pywcs.WCS(naxis=2)
        wcs.wcs.crpix = [0.5 * (nx + 1)] * 2
        wcs.wcs.crval = self.parameters["sky_center"]
        wcs.wcs.cdelt = [-dtheta, dtheta]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.cunit = ["deg"] * 2

        n_events = 0

        e = []
        x = []
        y = []
        for fn in self.filenames:
            with h5py.File(fn, "r") as f:
                d = f["data"]
                xx, yy = wcs.wcs_world2pix(d["xsky"][:], d["ysky"][:], 1)
                keepx = np.logical_and(xx >= 0.5, xx <= float(nx) + 0.5)
                keepy = np.logical_and(yy >= 0.5, yy <= float(nx) + 0.5)
                keep = np.logical_and(keepx, keepy)
                n_events += keep.sum()
                e.append(d["eobs"][keep])
                x.append(xx[keep])
                y.append(yy[keep])

        mylog.info(
            "Threw out %d events because they fell outside the " "field of view.",
            self.tot_num_events - n_events,
        )

        col_e = fits.Column(
            name="ENERGY", format="E", unit="eV", array=np.concatenate(e) * 1000.0
        )
        col_x = fits.Column(name="X", format="D", unit="pixel", array=np.concatenate(x))
        col_y = fits.Column(name="Y", format="D", unit="pixel", array=np.concatenate(y))

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
        tbhdu.header["TCRPX2"] = 0.5 * (nx + 1)
        tbhdu.header["TCRPX3"] = 0.5 * (nx + 1)
        tbhdu.header["TLMIN2"] = 0.5
        tbhdu.header["TLMIN3"] = 0.5
        tbhdu.header["TLMAX2"] = float(nx) + 0.5
        tbhdu.header["TLMAX3"] = float(nx) + 0.5
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

    def write_to_simput(self, prefix, overwrite=False):
        r"""
        Write events to a SIMPUT catalog that may be utilized by various
        instrument simulators.

        Parameters
        ----------
        prefix : string
            The filename prefix. The files to be written have the
            signature:
            f"{prefix}_simput.fits", f"{prefix}_phlist.fits", etc.
        overwrite : boolean, optional
            Set to True to overwrite previous files.
        """
        import unyt as u
        from astropy.coordinates import SkyCoord
        from soxs.simput import SimputCatalog, SimputPhotonList

        simput_file = f"{prefix}_simput.fits"

        begin_cat = True
        for i, fn in enumerate(self.filenames):
            if len(self.filenames) == 1:
                phlist_file = f"{prefix}_phlist.fits"
                name = os.path.basename(prefix)
            else:
                phlist_file = f"{prefix}_phlist.{i:04d}.fits"
                name = f"{os.path.basename(prefix)}.{i:04d}"
            with h5py.File(fn, "r") as f:
                d = f["data"]
                if d["eobs"].shape[0] > 0:
                    flux = (
                        np.sum(d["eobs"][()] * u.keV).to_value("erg")
                        / self.parameters["exp_time"]
                        / self.parameters["area"]
                    )

                    if self.observer == "internal":
                        c = SkyCoord(
                            d["xsky"][()], d["ysky"][()], unit="deg", frame="galactic"
                        )
                        ra = c.icrs.ra.value
                        dec = c.icrs.dec.value
                    else:
                        ra = d["xsky"][()]
                        dec = d["ysky"][()]
                    src = SimputPhotonList(ra, dec, d["eobs"][()], flux, name=name)

                    if begin_cat:
                        cat = SimputCatalog.from_source(
                            simput_file,
                            src,
                            src_filename=phlist_file,
                            overwrite=overwrite,
                        )
                        begin_cat = False
                    else:
                        cat.append(src, src_filename=phlist_file, overwrite=overwrite)
                else:
                    mylog.warning("No events found in file %s, so skipping.", fn)

    def write_fits_image(
        self, imagefile, fov, nx, emin=None, emax=None, overwrite=False
    ):
        r"""
        Generate an image by binning X-ray counts and write it to a FITS file.

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
            emin = -1.0
        if emax is None:
            emax = 1.0e10

        dtheta = fov.to_value("deg") / nx

        xbins = np.linspace(0.5, float(nx) + 0.5, nx + 1, endpoint=True)
        ybins = np.linspace(0.5, float(nx) + 0.5, nx + 1, endpoint=True)

        wcs = pywcs.WCS(naxis=2)
        wcs.wcs.crpix = [0.5 * (nx + 1)] * 2
        wcs.wcs.crval = self.parameters["sky_center"]
        wcs.wcs.cdelt = [-dtheta, dtheta]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.cunit = ["deg"] * 2

        H = np.zeros((nx, nx))

        for fn in self.filenames:
            with h5py.File(fn, "r") as f:
                d = f["data"]
                mask = np.logical_and(d["eobs"][:] >= emin, d["eobs"][:] <= emax)
                xx, yy = wcs.wcs_world2pix(d["xsky"][mask], d["ysky"][mask], 1)
                H += np.histogram2d(xx, yy, bins=[xbins, ybins])[0]

        hdu = fits.PrimaryHDU(H.T)

        hdu.header["MTYPE1"] = "EQPOS"
        hdu.header["MFORM1"] = "RA,DEC"
        hdu.header["CTYPE1"] = "RA---TAN"
        hdu.header["CTYPE2"] = "DEC--TAN"
        hdu.header["CRPIX1"] = 0.5 * (nx + 1)
        hdu.header["CRPIX2"] = 0.5 * (nx + 1)
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
        Bin observer-frame event energies into a spectrum and write it to
        a FITS binary table. This is for an *unconvolved* spectrum.

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
        spec = np.zeros(nchan)
        ebins = np.linspace(emin, emax, nchan + 1, endpoint=True)
        emid = 0.5 * (ebins[1:] + ebins[:-1])

        for fn in self.filenames:
            with h5py.File(fn, "r") as f:
                d = f["data"]
                spec += np.histogram(d["eobs"][:], bins=ebins)[0]

        col1 = fits.Column(
            name="CHANNEL", format="1J", array=np.arange(nchan).astype("int32") + 1
        )
        col2 = fits.Column(name="ENERGY", format="1D", array=emid.astype("float64"))
        col3 = fits.Column(name="COUNTS", format="1J", array=spec.astype("int32"))
        col4 = fits.Column(
            name="COUNT_RATE", format="1D", array=spec / self.parameters["exp_time"]
        )

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
