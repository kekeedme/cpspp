"""This file contains applications of the integral functions written in the
integrals.py file"""
#Don't necessarily run all the programs if you don't need them, it will make the program slow
import numpy as np
from src.useful_math_functions.integrals import function_integrationtz
from src.useful_math_functions.integrals import adaptivetz
from src.useful_math_functions.integrals import Rombergtz
from src.useful_math_functions.integrals import gaussquad
from src.useful_math_functions.integrals import function_integrationsimp
from pylab import plot, show, imshow, xlim, legend, title, hot, show

#==============================================================================================
# Application of integration to compute Bessel functions

# We will need to calculate the Bessel functions (integrals of trig functions)
# We will use the simpson integrator

npoints = 1000
x_vals = np.linspace(0, 20, npoints + 1)

# Math function J_m(x)=1/pi \int{cos(m\theta-xsin\theta)d\theta}


# defining the inner function which needs to integrated. Note that the values
# m and xval will be taken from the loop, and thus does not need to be specified
# when we call the function
def bessel_trig(theta):
    return np.cos(m_val * theta - xval * np.sin(theta))


# the loops which will carry out the integrations over values of m and x
for m_val in range(0, 3):
    j_vals = (
        []
    )  # initialize list which will
    # contain values of bessel func and reinitialize after each loop over m
    for xval in x_vals:
        # performing integration over values of x from 0 to 20 for fixed m
        jval = (1 / np.pi) * function_integrationsimp(bessel_trig, 0, np.pi, npoints)

        j_vals.append(jval)
    plot(x_vals, j_vals, label="m = " + str(m_val))
xlim(x_vals[0], x_vals[-1])
title("Bessel function")
legend()
show()

# ==============================================================================================
# Diffraction limit exercise (uses functions from previous exercise)
# We will make a density plot of the intensity of the light in the diffraction pattern
# from the Bessel function calculated above, we will calculate the intensity
# in order to generate the diffraction pattern

# math function = I(r)=(J_1(kr)**2/(kr))

# we need to calculate the bessel function for m = 1 and xval = kr
# k =2pi/\lambda

wavelength = 5 * 10 ** (-5)  # cm \lambda
k_bessel = np.pi / wavelength  # k in formula
m_val = 1  # fixing m equal to 1
points = 200  # number of grid points along each side
bound = 2 ** (-0.5) * 10 ** (-4)
xpoints = np.linspace(-bound, bound, points)  # array of radius values in cm
ypoints = xpoints
x_mesh, y_mesh = np.meshgrid(xpoints, ypoints)
r_vals = np.sqrt(x_mesh**2 + y_mesh**2)
kr_vals = k_bessel * r_vals  # Calculating the argument in the function
intensity = np.empty([points, points])  # making an array to store intensity
for i in range(points):
    for j in range(points):
        xval = kr_vals[i, j]
        # getting the integration values for each x values at different grid points
        jval = (1 / np.pi) * function_integrationsimp(bessel_trig, 0, np.pi, points)

        intensity[i, j] = (jval**2) / (xval)

imshow(intensity, vmax=0.05, origin="lower", extent=[-bound, bound, -bound, bound])
hot()
show()

# ==============================================================================================
