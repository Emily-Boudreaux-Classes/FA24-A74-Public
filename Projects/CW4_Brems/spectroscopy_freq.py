'''
This module contains some useful tools for spectroscopy.
Right now, that's basically just the Planck function.
In general, it can be *very* helpful to store commonly
used code tools in modules (*.py files) that you can
import into notebooks, or other modules.
'''

import numpy as np
from astropy import units, constants
import astropy.io.ascii
import scipy.interpolate

def planck_intensity(frequency, temperature):
    '''
    This function calculates the thermal emission intensity spectrum of a surface.

        Inputs:
            wavelength = numpy array of wavelengths (with astropy units)
            temperature = a single number, the temperature (with astropy units)

        Outputs:
            Returns an array of thermal emission intensities,
            in astropy units of ergs/(cm^2*Hz*sr). This is a flux, which has
            already been integrated over solid angle.
    '''

    # define variables as shortcut to the constants we need
    h = constants.h
    k = constants.k_B
    c = constants.c

    # this is the thing that goes into the exponent (it's units better cancel!)
    u = (h*frequency/(k*temperature)).decompose()
    try:
        print(u[0])
    except:
        pass

    # calculate the intensity from the Planck function
    intensity = (2*h*frequency**3/c**2/(np.exp(u) - 1))/units.steradian
    try:
        print(intensity[0])
    except:
        pass

    # return the intensity
    return intensity.to('erg/(cm**2*s*Hz*sr)')


def planck_flux(frequency, temperature):
    '''
    This function calculates the thermal emission flux spectrum of a surface.

        Inputs:
            wavelength = numpy array of wavelengths (with astropy units)
            temperature = a single number, the temperature (with astropy units)

        Outputs:
            Returns an array of thermal emission fluxes,
            in astropy units of W/(m^2*micron). This is a flux, which has
            already been integrated over solid angle.
    '''

    # calculate the flux, knowing the angle integral will be pi steradians (for isotropic emission)
    flux = planck_intensity(frequency, temperature)*np.pi*units.steradian

    # return the flux, in convenient units
    return flux.to('W/(m**2*Hz)')
