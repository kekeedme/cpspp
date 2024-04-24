"""This file will also contain instances where we apply our integration function to perform
relevant physics calculations"""

import numpy as np
from pylab import plot, show, xlabel, ylabel, xlim, ylim, title
from src.useful_math_functions.integrals import gaussquad


# ================================================================================================
# PERIOD of an anharmonic oscillator
# taken from exercise 5.10 (p173) from Newman Computational Physics


# We will compute the period of an anharmonic oscillator for a potential of V(x) = x^4
# a t= 0 x = a (max amplitude) and when x = 0, T= (1/4)*T the first time
# def invapotential(xval):
#    """This function calculate the integrand
#    math:1/sqrt(V(a)-V(x))
#    :param: xval only the values of x integration"""
#    return np.sqrt(8 * MASS) * 1 / (np.sqrt((amps**4) - (xval**4)))

# ================================================================================================
# DIFFRACTION OF PLANE WAVE #Taken from exercise 5.11 (p174) from Newman Computational Physics


def ufunc(xval):
    """This function is related to the wavelength and the position z
    :param: xval the height value"""
    return xval * np.sqrt(2 / (WAVELENGT * ZVAL))


def cufunc(time):
    """This function is the integrand of the C(U) function it will be a cosine function
    :param: time"""
    return np.cos(0.5 * np.pi * time**2)


def sufunc(time):
    """This function is the integrand of the S(U) function it will be a sine function
    :param: time"""
    return np.sin(0.5 * np.pi * time**2)


# ================================================================================================

# Below starts the if __name__ = "__main__" to use the functions written above

if __name__ == "__main__":
    #### for the anharmonic oscillator exercise###
    # MASS = 1  # defining the mass constant
    # we will integrate for a = 0, to a = 2, and make a plot
    # period_list = []  # list that will contain the values of the calculated periods
    # amplitude_list = np.linspace(0.1, 2, 50)  # list that contains the amplitudes

    # for amplitude in amplitude_list:
    #    amps = amplitude
    #    period_list.append(
    #        gaussquad(invapotential, 0, amplitude, 20)
    #    )  # performing the integral
    # and appending the values to the list
    # plot(amplitude_list, period_list)
    # xlim(amplitude_list[0], amplitude_list[-1])
    # ylim(period_list[-1], period_list[0])
    # title("Period of anharmonic oscillator")
    # xlabel("Amplitude")
    # ylabel("Period")
    # show()

    #### For the sound diffraction exercise###
    WAVEL = 1  # in meters
    WAVELENGT = WAVEL
    _Z = 3  # in meters
    ZVAL = _Z

    # Recall math: I/I_0 = (1/8) * ([2C(u)+1]^2 + [2S(u)+1]^2)

    xlist = np.linspace(-5, 5, 100)  # list of xvals
    iratio = []  # list of intensity ratio I/I_0
    # writing for loop to calculate I/I_0 for a range of x -5 to 5 m

    for xvals in xlist:
        ufuncval = ufunc(xvals)  # calculate the ufunction at the given xval
        cuvals = gaussquad(
            cufunc, 0, ufuncval, 50
        )  # integrate from 0 to ufuncval using gaussquad
        suvals = gaussquad(
            sufunc, 0, ufuncval, 50
        )  # integrate from 0 to ufuncval using gaussquad
        iratioval = (1 / 8) * ((2 * cuvals + 1) ** 2 + (2 * suvals + 1) ** 2)
        iratio.append(iratioval)

    # making the plot of iratio as a function of x
    plot(xlist, iratio)
    xlim(-5, 5)
    ylim(iratio[0], 1.4 * iratio[-1])
    xlabel("X positions")
    ylabel("Ratio of intensity")
    title("Ratio of intensity versus position")
    show()
