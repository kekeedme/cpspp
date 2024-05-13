"""This file contains one application of a numerical two-dimensional integral using
gaussian quadrature"""

import numpy as np
from pylab import show, xlabel, ylabel, title, plot, xlim, legend
from scipy.constants import gravitational_constant as _G

from src.useful_math_functions.integrals import gaussquad2d
from src.useful_math_functions.integrals import gaussquad2_2d

# defining some of the constants we will use
SIGMA = 1000 / 100  # mass per unit area total mass in kg divided by total area in m^2
MASS = 1  # mass of the object (in Kg) experiencing gravity from the sheet

# Lists that are needed
z_list = np.linspace(0, 10, 100)  # array of distances of the object from the sheet


def dforce_z(xval, yval, zval):
    """This function is the infinitesimal gravitational force along the z axis
    experienced by an object of mass m, due to a sheet of mass density sigma
    :param: xval yval,zval the values of the x,y,z coordinates"""
    return (_G * SIGMA * zval * MASS) / (xval**2 + yval**2 + zval**2) ** (3 / 2)


# We will perform two changes of variables to get a function which we can integrate to z->0


def polarforce(theta, phi):
    """This function is the infinitesimal gravitational force along the z axis when z is zero
    experienced by an object of mass m, due to a sheet of mass density sigma
    :param: theta and phi which are the variable used in the sub to make the intg finite
    """
    deno = (1 / (np.cos(theta) ** 2 * np.cos(phi) ** 2)) * (
        1 / (np.tan(theta) ** 2 + np.tan(phi) ** 2 + 1) ** (3 / 2)
    )
    num = _G * SIGMA * MASS
    return num * deno


def exactforce(zval):
    """This function calculates the force from a finite disk or (sheet of the same mass density)
    as a function of the z-distance from the sheet
    note: here we took the radius of the disk such that the area is equal to the area of the
    finite sheet so rˆ2 = Lˆ2/pi = 31.4"""
    return (2 * _G * 1000 * MASS / 31.4) * ((-zval / np.sqrt(zval**2 + 31.4)) + 1)


force_list = []  # making an empty list that will contain the force values
force = gaussquad2_2d(polarforce, -np.pi / 2, np.pi / 2, 100)
force_list.append(
    force
)  # calculate the first point using integral on infinite range since z =0

# Calculating for other points
for zvalues in z_list[1:]:
    force = gaussquad2d(dforce_z, -5, 5, zvalues, 200)
    force_list.append(force)

# calculating the exact values
exact = exactforce(z_list)

# making the plot
plot(z_list, force_list, ".", label="Numerical")
plot(z_list, exact, label="Exact")
legend()
xlabel("Position from the center of the sheet (m)")
xlim(z_list[0] - 0.1, z_list[-1] + 0.1)
ylabel("Force (N)")
title("Gravitational force due to uniform sheet as a function of altitude")
show()
