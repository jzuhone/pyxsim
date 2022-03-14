import numpy as np
from yt.utilities.cosmology import Cosmology


def make_xray_fields(ds, source_model, redshift=0.0,
                     dist=None, cosmology=None):

    spectral_norm = 1.0

    source_model.setup_model(ds, redshift, spectral_norm)

    ftype = source_model.ftype

    def _emissivity_field(field, data):
        ret = data.ds.arr(source_model.process_data("energy_field", data),
                          "keV/s")
        return ret*data[ftype, "density"]/data[ftype, "mass"]

    if hasattr(source_model, "emin"):
        emiss_name = (
            ftype, f"xray_emissivity_{source_model.emin.value}_{source_model.emax.value}_keV"
        )
        emiss_dname = rf"\epsilon_{{X}} ({source_model.emin.value}-{source_model.emax.value} keV)"
    else:
        emiss_name = (
            ftype, f"xray_emissivity_{source_model.e0.value}_keV"
        )
        emiss_dname = rf"\epsilon_{{X}} ({source_model.e0.value} keV)"

    ds.add_field(
        emiss_name,
        function=_emissivity_field,
        display_name=emiss_dname,
        sampling_type="local",
        units="erg/cm**3/s",
    )

    def _luminosity_field(field, data):
        return data.ds.arr(source_model.process_data("energy_field", data),
                           "keV/s")

    lum_name = (
        ftype, emiss_name[1].replace("emissivity", "luminosity")
    )
    lum_dname = emiss_dname.replace("\epsilon", "\rm{{L}}") 

    ds.add_field(
        lum_name,
        function=_luminosity_field,
        display_name=lum_dname,
        sampling_type="local",
        units="erg/s",
    )

    def _photon_emissivity_field(field, data):
        ret = data.ds.arr(source_model.process_data("photon_field", data),
                          "photons/s")
        return ret * data[ftype, "density"] / data[ftype, "mass"]

    phot_emiss_name = (
        ftype, emiss_name[1].replace("emissivity", "photon_emissivity")
    )

    ds.add_field(
        phot_emiss_name,
        function=_photon_emissivity_field,
        display_name=emiss_dname,
        sampling_type="local",
        units="photons/cm**3/s",
    )

    xray_fields = [emiss_name, lum_name, phot_emiss_name]

    if redshift > 0.0 or dist is not None:

        if dist is None:
            if cosmology is None:
                if hasattr(ds, "cosmology"):
                    cosmology = ds.cosmology
                else:
                    cosmology = Cosmology()
            D_L = cosmology.luminosity_distance(0.0, redshift)
            angular_scale = 1.0 / cosmology.angular_scale(0.0, redshift)
            dist_fac = ds.quan(
                1.0 / (4.0 * np.pi * D_L * D_L * angular_scale * angular_scale).v,
                "rad**-2",
                )
        else:
            redshift = 0.0  # Only for local sources!
            try:
                # normal behaviour, if dist is a YTQuantity
                dist = ds.quan(dist.value, dist.units)
            except AttributeError as e:
                try:
                    dist = ds.quan(*dist)
                except (RuntimeError, TypeError):
                    raise TypeError(
                        "dist should be a YTQuantity or a (value, unit) tuple!"
                    ) from e

            angular_scale = dist / ds.quan(1.0, "radian")
            dist_fac = ds.quan(
                1.0 / (4.0 * np.pi * dist * dist * angular_scale * angular_scale).v,
                "rad**-2",
                )

        ei_name = (
            ftype, emiss_name[1].replace("emissivity", "intensity")
        )
        ei_dname = emiss_name.replace("\epsilon", "I")
        def _intensity_field(field, data):
            I = dist_fac * data[emiss_name]
            return I.in_units("erg/cm**3/s/arcsec**2")

        ds.add_field(
            ei_name,
            function=_intensity_field,
            display_name=ei_dname,
            sampling_type="local",
            units="erg/cm**3/s/arcsec**2",
        )

        i_name = (
            ftype, phot_emiss_name[1].replace("emissivity", "intensity")
        )
        i_dname = phot_emiss_name.replace("\epsilon", "I")

        def _photon_intensity_field(field, data):
            I = (1.0 + redshift) * dist_fac * data[phot_emiss_name]
            return I.in_units("photons/cm**3/s/arcsec**2")

        ds.add_field(
            i_name,
            function=_photon_intensity_field,
            display_name=i_dname,
            sampling_type="local",
            units="photons/cm**3/s/arcsec**2",
        )

        xray_fields += [ei_name, i_name]

    return xray_fields