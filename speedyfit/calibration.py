import os

from astropy.io import fits

from speedyfit.filters import synthetic_flux

basedir = os.path.dirname(__file__)

alpha_lyr_path = 'calibration/VEGA_alpha_lyr_stis_008.fits'


def get_vega_calibrator():
    """
    Returns the wavelength and flux of Vega / Alpha Lyr
    Wavelength in AA
    Flux in erg/s/cm2/AA
    """
    filepath = os.path.join(basedir, alpha_lyr_path)

    hdu = fits.open(filepath)
    data = hdu[1].data
    hdu.close()

    wave, flux = data['WAVELENGTH'], data['FLUX']

    return wave, flux


def get_Flambda(photband):
    """
    Get F0 (Flambda) for a given photometric pass band.

    Returns F0 in erg/s/cm2/AA

    Vegamag = -2.5 Log F_lam - C with F_lam in erg/s/cm2/AA
    Flux computed as 10**(-meas/2.5)*F0
    Magnitude computed as -2.5*log10(Fmeas/F0)
    """
    # -- get calibrator
    wave, flux = get_vega_calibrator()

    # -- calculate synthetic fluxes
    syn_flux = synthetic_flux(wave, flux, [photband])

    return syn_flux

# def calibrate_photband():
#     """
#     Calibrate photometry.
#
#     Vegamag = -2.5 Log F_lam - C with F_lam in erg/s/cm2/AA
#     Flux computed as 10**(-meas/2.5)*F0
#     Magnitude computed as -2.5*log10(Fmeas/F0)
#     """
#
#     # -- get calibrator
#     wave, flux = get_vega_calibrator()
#
#     # -- calculate synthetic fluxes
#     syn_flux = synthetic_flux(wave, flux, zp['photband'])
#
#     # -- set Flam0 for all the bands
#     zp['Flam0'] = syn_flux
#     zp['Flam0_units'] = 'erg/s/cm2/AA'
#
#     # -- set the central wavelengths of the bands
#     zp['eff_wave'] = filters.eff_wave(zp['photband'])
#
#     return zp