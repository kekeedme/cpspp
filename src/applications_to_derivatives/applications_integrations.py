"""This file contains applications of the integral functions written in the
integrals.py file"""

# Don't necessarily run all the programs if you don't need them, it will make the program slow

import numpy as np

# from src.useful_math_functions.integrals import function_integrationtz
# from src.useful_math_functions.integrals import adaptivetz
# from src.useful_math_functions.integrals import rombergtz
from src.useful_math_functions.integrals import gaussquad

# from src.useful_math_functions.integrals import function_integrationsimp
from pylab import plot, xlim, ylim, xlabel, ylabel, title, show

# hot, imshow,legend

# 1=============================================================================================
# Application of integration to compute Bessel functions

# We will need to calculate the Bessel functions (integrals of trig functions)
# We will use the simpson integrator

# npoints = 1000
# x_vals = np.linspace(0, 20, npoints + 1)

# Math function J_m(x)=1/pi \int{cos(m\theta-xsin\theta)d\theta}


# defining the inner function which needs to integrated. Note that the values
# m and xval will be taken from the loop, and thus does not need to be specified
# when we call the function
# def bessel_trig(theta):
#    return np.cos(m_val * theta - xval * np.sin(theta))


# the loops which will carry out the integrations over values of m and x
# for m_val in range(0, 3):
#   j_vals = (
#      []
# )  # initialize list which will
# contain values of bessel func and reinitialize after each loop over m
# for xval in x_vals:
# performing integration over values of x from 0 to 20 for fixed m
#   jval = (1 / np.pi) * function_integrationsimp(bessel_trig, 0, np.pi, npoints)

#  j_vals.append(jval)
# plot(x_vals, j_vals, label="m = " + str(m_val))
# xlim(x_vals[0], x_vals[-1])
# title("Bessel function")
# legend()
# show()

# 2 =============================================================================================
# DIFFRACTION limit exercise (uses functions from previous exercise)
# We will make a density plot of the intensity of the light in the diffraction pattern
# from the Bessel function calculated above, we will calculate the intensity
# in order to generate the diffraction pattern

# math function = I(r)=(J_1(kr)**2/(kr))

# we need to calculate the bessel function for m = 1 and xval = kr
# k =2pi/\lambda

# wavelength = 5 * 10 ** (-5)  # cm \lambda
# k_bessel = np.pi / wavelength  # k in formula
# m_val = 1  # fixing m equal to 1
# points = 200  # number of grid points along each side
# bound = 2 ** (-0.5) * 10 ** (-4)
# xpoints = np.linspace(-bound, bound, points)  # array of radius values in cm
# ypoints = xpoints
# x_mesh, y_mesh = np.meshgrid(xpoints, ypoints)
# r_vals = np.sqrt(x_mesh**2 + y_mesh**2)
# kr_vals = k_bessel * r_vals  # Calculating the argument in the function
# intensity = np.empty([points, points])  # making an array to store intensity
# for i in range(points):
#   for j in range(points):
#      xval = kr_vals[i, j]
# getting the integration values for each x values at different grid points
#     jval = (1 / np.pi) * function_integrationsimp(bessel_trig, 0, np.pi, points)

#    intensity[i, j] = (jval**2) / (xval)

# imshow(intensity, vmax=0.05, origin="lower", extent=[-bound, bound, -bound, bound])
# hot()
# show()

# 3==============================================================================================
# HEAT CAPACITY OF SOLIDS

# recall math: C_v = 9V\rho k_B (T/theta)^3 *
#   integral (x^4e^x)/(e^x-1)^2 (eval from 0 to upperbound), where x = upperbound (theta/T)


# writing the polynomial which will need to be integrated to calc the heat capacity
def polynomial(temp):
    """This function will be integrated from 0 to the upperbound specified
    in the heat capacity function
    :param: temp which is the Debye temperature THETA divided by the temperature of interests
    """
    return ((temp**4) * np.exp(temp)) / ((np.exp(temp) - 1) ** 2)


def heatcapacity(temperature):
    """This function will calculate the heat capacity for a given solid
    using the Debye formula, where the integration will be performed using
    gauss quadrature
    :param: temperature, the temperature at which we wish to calculate the heat capacity
    """
    upperbound = THETA / temperature  # defining the upperbound of the integral
    integral_result = gaussquad(
        polynomial, 0, upperbound, 50
    )  # performing integration for N=50
    heatcapacityval = (
        9 * RHO * VOLUME * KBOLTZMAN * ((temperature / THETA) ** 3) * integral_result
    )
    return heatcapacityval


if __name__ == "__main__":
    # defining the necessary variables
    RHO = 6.022e28  # number density of atoms
    THETA = 428  # Debye temperature in Kelvin
    KBOLTZMAN = 1.380649e-23
    VOLUME = 0.001  # in cubic meters

    temp_array = np.linspace(5, 500, 100)
    heat_capa = []
    for k, temperat in enumerate(temp_array):
        heat_capacity = heatcapacity(temperat)
        heat_capa.append(heat_capacity)
    plot(temp_array, heat_capa)
    xlim(temp_array[0], temp_array[-1])
    ylim(heat_capa[0], heat_capa[-1])
    xlabel("Temperature values in K")
    ylabel("Heat Capacity in J/K ")
    title("Heat Capacity of Aluminum")
    show()
