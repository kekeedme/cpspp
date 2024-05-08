"""This file contains one application of a numerical two-dimensional integral using
gaussian quadrature"""

import numpy as np
from pylab import show, xlabel, ylabel, title, plot
from scipy.constants import gravitational_constant as _G

from src.useful_math_functions.integrals import gaussquad2d

# defining some of the constants we will use
SIGMA = 1000 / 100  # mass per unit area total mass in kg divided by total area in m^2
MASS = 1  # mass of the object (in Kg) experiencing gravity from the sheet

#Lists that are needed
z_list = np.linspace(0, 10, 100)  # array of distances of the object from the sheet
force_list = []


def dforce_z(xval, yval, zval):
    """This function is the infinitesimal gravitational force along the z axis
    experienced by an object of mass m, due to a sheet of mass density sigma
    :param: xval yval,zval the values of the x,y,z coordinates"""
    return (_G * SIGMA * zval * MASS) / (xval ** 2 + yval ** 2 + zval ** 2) ** (3 / 2)


for zvalues in z_list:
    force = gaussquad2d(dforce_z, -5, 5, zvalues, 100)
    force_list.append(force)

# making the plot
plot(z_list, force_list, "--")
xlabel("Position from the center of the sheet (m)")
ylabel("Force (N)")
title("Gravitational force due to uniform sheet as a function of altitude")
show()
