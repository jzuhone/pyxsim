"""
Classes for generating lists of detected events
"""
import numpy as np
from pyxsim.utils import mylog
try:
    from yt.visualization.fits_image import assert_same_wcs
except ImportError:
    from yt.utilities.fits_image import assert_same_wcs
from yt.units.yt_array import YTQuantity, YTArray, uconcatenate
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import h5py
from pyxsim.utils import force_unicode, validate_parameters, parse_value, \
    ParameterDict
from soxs.simput import write_photon_list
from soxs.instrument import RedistributionMatrixFile
import os
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    communication_system, parallel_capable, get_mpi_type

old_parameter_keys = {"ExposureTime": "exp_time",
                      "Area": "area",
                      "Redshift": "redshift",
                      "AngularDiameterDistance": "d_a",
                      "RMF": "rmf",
                      "ARF": "arf",
                      "Telescope": "telescope",
                      "Instrument": "instrument",
                      "Mission": "mission",
                      "ChannelType": "channel_type"}

comm = communication_system.communicators[-1]

def communicate_events(my_events, root=0):
    if parallel_capable:
        new_events = {}
        mpi_int = get_mpi_type("int32")
        mpi_double = get_mpi_type("float64")
        local_num_events = my_events["xpix"].size
        sizes = comm.comm.gather(local_num_events, root=root)
        if comm.rank == 0:
            num_events = sum(sizes)
            disps = [sum(sizes[:i]) for i in range(len(sizes))]
            for key in my_events:
                if key in ["pi", "pha"]:
                    dtype = "int32"
                else:
                    dtype = "float64"
                new_events[key] = np.zeros(num_events, dtype=dtype)
        else:
            sizes = []
            disps = []
            for key in my_events:
                new_events[key] = np.empty([])
        for key in my_events:
            if key in ["pi", "pha"]:
                mpi_type = mpi_int
            else:
                mpi_type = mpi_double
            comm.comm.Gatherv([my_events[key], local_num_events, mpi_type],
                              [new_events[key], (sizes, disps), mpi_type], root=root)
            if key == "eobs":
                new_events[key] = YTArray(new_events[key], "keV")
            if key.endswith("sky"):
                new_events[key] = YTArray(new_events[key], "deg")
        return new_events
    else:
        return my_events

class EventList(object):

    def __init__(self, events, parameters, wcs=None):
        self.events = events
        self.parameters = ParameterDict(parameters, "EventList", old_parameter_keys)
        self.num_events = comm.mpi_allreduce(events["xpix"].shape[0])
        if wcs is None:
            self.wcs = pywcs.WCS(naxis=2)
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
        return type(self)(events, self.parameters, wcs=self.wcs)

    def __iter__(self):
        return iter(self.events)

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
        return type(self)(new_events, self.parameters, wcs=self.wcs)

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

        num_events = d["xpix"].size
        start_e = comm.rank*num_events//comm.size
        end_e = (comm.rank+1)*num_events//comm.size

        events["xpix"] = d["xpix"][start_e:end_e]
        events["ypix"] = d["ypix"][start_e:end_e]
        events["eobs"] = YTArray(d["eobs"][start_e:end_e], "keV")
        if "rmf" in p:
            parameters["rmf"] = force_unicode(p["rmf"].value)
            parameters["arf"] = force_unicode(p["arf"].value)
            parameters["channel_type"] = force_unicode(p["channel_type"].value)
            parameters["mission"] = force_unicode(p["mission"].value)
            parameters["telescope"] = force_unicode(p["telescope"].value)
            parameters["instrument"] = force_unicode(p["instrument"].value)
            chantype = parameters["channel_type"]
            events[chantype] = d[chantype][start_e:end_e]

        f.close()

        if "rmf" in p:
            return ConvolvedEventList(events, parameters)
        else:
            return EventList(events, parameters)

    @classmethod
    def from_fits_file(cls, fitsfile):
        """
        Initialize an :class:`~pyxsim.event_list.EventList` from a FITS 
        file with filename *fitsfile*.
        """
        hdulist = pyfits.open(fitsfile, memmap=True)

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

        num_events = tblhdu.header["NAXIS2"]
        start_e = comm.rank*num_events//comm.size
        end_e = (comm.rank+1)*num_events//comm.size

        events["xpix"] = tblhdu.data["X"][start_e:end_e]
        events["ypix"] = tblhdu.data["Y"][start_e:end_e]
        events["eobs"] = YTArray(tblhdu.data["ENERGY"][start_e:end_e]/1000., "keV")

        if "rmf" in tblhdu.header:
            chan_type = parameters["channel_type"].upper()
            parameters["rmf"] = tblhdu.header["RMF"]
            parameters["arf"] = tblhdu.header["ARF"]
            parameters["channel_type"] = tblhdu.header["CHANTYPE"]
            parameters["mission"] = tblhdu.header["MISSION"]
            parameters["telescope"] = tblhdu.header["TELESCOP"]
            parameters["instrument"] = tblhdu.header["INSTRUME"]
            events[chan_type] = tblhdu.data[chan_type][start_e:end_e]

        hdulist.close()

        if "rmf" in tblhdu.header:
            return ConvolvedEventList(events, parameters)
        else:
            return EventList(events, parameters)

    def write_fits_file(self, fitsfile, overwrite=False):
        """
        Write events to a FITS binary table file with filename *fitsfile*.
        Set *overwrite* to True if you need to overwrite a previous file.
        """
        from astropy.time import Time, TimeDelta

        events = communicate_events(self.events)

        if comm.rank == 0:

            exp_time = float(self.parameters["exp_time"])

            t_begin = Time.now()
            dt = TimeDelta(exp_time, format='sec')
            t_end = t_begin + dt

            col_e = pyfits.Column(name='ENERGY', format='E', unit='eV',
                                  array=events["eobs"].in_units("eV").d)
            col_x = pyfits.Column(name='X', format='D', unit='pixel',
                                  array=events["xpix"])
            col_y = pyfits.Column(name='Y', format='D', unit='pixel',
                                  array=events["ypix"])

            cols = [col_e, col_x, col_y]

            if "channel_type" in self.parameters:
                chantype = self.parameters["channel_type"]
                if chantype == "pha":
                    cunit = "adu"
                elif chantype == "pi":
                    cunit = "Chan"
                col_ch = pyfits.Column(name=chantype.upper(), format='1J',
                                       unit=cunit, array=events[chantype])
                cols.append(col_ch)

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
            if "channel_type" in self.parameters:
                rmf = RedistributionMatrixFile(self.parameters["rmf"])
                tbhdu.header["TLMIN4"] = rmf.cmin
                tbhdu.header["TLMAX4"] = rmf.cmax
                tbhdu.header["RESPFILE"] = os.path.split(self.parameters["rmf"])[-1]
                tbhdu.header["PHA_BINS"] = rmf.n_ch
                tbhdu.header["ANCRFILE"] = os.path.split(self.parameters["arf"])[-1]
                tbhdu.header["CHANTYPE"] = self.parameters["channel_type"]
                tbhdu.header["MISSION"] = self.parameters["mission"]
                tbhdu.header["TELESCOP"] = self.parameters["telescope"]
                tbhdu.header["INSTRUME"] = self.parameters["instrument"]
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
            tbhdu.header["DATE"] = t_begin.tt.isot
            tbhdu.header["DATE-OBS"] = t_begin.tt.isot
            tbhdu.header["DATE-END"] = t_end.tt.isot

            hdulist = [pyfits.PrimaryHDU(), tbhdu]

            if "channel_type" in self.parameters:
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

            pyfits.HDUList(hdulist).writeto(fitsfile, overwrite=overwrite)

        comm.barrier()

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
        events = communicate_events(self.events)

        if comm.rank == 0:

            if emin is None:
                emin = events["eobs"].min().value
            if emax is None:
                emax = events["eobs"].max().value

            idxs = np.logical_and(events["eobs"].d >= emin, events["eobs"].d <= emax)
            flux = np.sum(events["eobs"][idxs].in_units("erg")) / \
                   self.parameters["exp_time"]/self.parameters["area"]

            write_photon_list(prefix, prefix, flux.v, events["xsky"][idxs].d,
                              events["ysky"][idxs].d, events["eobs"][idxs].d,
                              overwrite=overwrite)

        comm.barrier()

    def write_h5_file(self, h5file):
        """
        Write an :class:`~pyxsim.event_list.EventList` to the HDF5 file given by *h5file*.
        """
        events = communicate_events(self.events)

        if comm.rank == 0:

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
            d.create_dataset("xpix", data=events["xpix"])
            d.create_dataset("ypix", data=events["ypix"])
            d.create_dataset("xsky", data=events["xsky"].d)
            d.create_dataset("ysky", data=events["ysky"].d)
            d.create_dataset("eobs", data=events["eobs"].d)
            if "rmf" in self.parameters:
                p.create_dataset("arf", data=self.parameters["arf"])
                p.create_dataset("rmf", data=self.parameters["rmf"])
                p.create_dataset("channel_type", data=self.parameters["channel_type"])
                p.create_dataset("mission", data=self.parameters["mission"])
                p.create_dataset("telescope", data=self.parameters["telescope"])
                p.create_dataset("instrument", data=self.parameters["instrument"])
                chantype = self.parameters["channel_type"]
                d.create_dataset(chantype, data=events[chantype])

            f.close()

        comm.barrier()

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

        if parallel_capable:
            H = comm.comm.reduce(H, root=0)

        if comm.rank == 0:

            hdu = pyfits.PrimaryHDU(H.T)

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

        comm.barrier()

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
        espec = self["eobs"].d
        spec, ee = np.histogram(espec, bins=nchan, range=(emin, emax))
        bins = 0.5*(ee[1:]+ee[:-1])

        if parallel_capable:
            spec = comm.comm.reduce(spec, root=0)

        if comm.rank == 0:

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

        comm.barrier()

class ConvolvedEventList(EventList):

    def write_simput_file(self, prefix, overwrite=False, emin=None, emax=None):
        mylog.error("Writing SIMPUT files is only supported if you didn't convolve with responses!")
        raise TypeError("Writing SIMPUT files is only supported if you didn't convolve with responses!")

    def write_channel_spectrum(self, specfile, overwrite=False):
        r"""
        Bin event channels into a spectrum and write it to a FITS binary table. 

        Parameters
        ----------
        specfile : string
            The name of the FITS file to be written.
        overwrite : boolean, optional
            Set to True to overwrite a previous file.
        """
        spectype = self.parameters["channel_type"]
        rmf = RedistributionMatrixFile(self.parameters["rmf"])
        minlength = rmf.n_ch
        if rmf.cmin == 1: 
            minlength += 1
        spec = np.bincount(self[spectype], minlength=minlength)
        if rmf.cmin == 1: 
            spec = spec[1:]
        bins = (np.arange(rmf.n_ch)+rmf.cmin).astype("int32")

        if parallel_capable:
            spec = comm.comm.reduce(spec, root=0)

        if comm.rank == 0:

            col1 = pyfits.Column(name='CHANNEL', format='1J', array=bins)
            col2 = pyfits.Column(name=spectype.upper(), format='1D', array=bins.astype("float64"))
            col3 = pyfits.Column(name='COUNTS', format='1J', array=spec.astype("int32"))
            col4 = pyfits.Column(name='COUNT_RATE', format='1D', array=spec/float(self.parameters["exp_time"]))

            coldefs = pyfits.ColDefs([col1, col2, col3, col4])

            tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
            tbhdu.name = "SPECTRUM"

            tbhdu.header["DETCHANS"] = spec.shape[0]
            tbhdu.header["TOTCTS"] = spec.sum()
            tbhdu.header["EXPOSURE"] = float(self.parameters["exp_time"])
            tbhdu.header["LIVETIME"] = float(self.parameters["exp_time"])
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
            tbhdu.header["RESPFILE"] = os.path.split(self.parameters["rmf"])[-1]
            tbhdu.header["ANCRFILE"] = os.path.split(self.parameters["arf"])[-1]
            tbhdu.header["MISSION"] = self.parameters["mission"]
            tbhdu.header["TELESCOP"] = self.parameters["telescope"]
            tbhdu.header["INSTRUME"] = self.parameters["instrument"]
            tbhdu.header["AREASCAL"] = 1.0
            tbhdu.header["CORRSCAL"] = 0.0
            tbhdu.header["BACKSCAL"] = 1.0

            hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])

            hdulist.writeto(specfile, overwrite=overwrite)

        comm.barrier()