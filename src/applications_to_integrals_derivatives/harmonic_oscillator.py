"""This file contains applications of integration to the qm harmonic oscillator"""

from math import factorial as fact
import numpy as np
from pylab import plot, show, title, xlabel, ylabel, legend

from src.useful_math_functions.integrals import gaussquad

# Harmonic oscillator
# taken from exercise 5.13 (p182) from Newman Computational Physics


def hermite(nval, xval):
    """This function will calculate the nth hermite polynomial for given xval
    using recursion calls
    :param: nval index value positive integer
    :param: xval values at which the polynomial should be evaluated"""
    if n < 0:
        print("The value of n should be greater than or equal to zero")
        return None
    if nval == 0:
        return 1
    if nval == 1:
        return 2 * xval
    # else
    return (2 * xval * hermite(nval - 1, xval)) - (
        2 * (nval - 1) * hermite(nval - 2, xval)
    )


def harmoscillator(nval, xval):
    """This function computes the nth quantum harmonic oscillator function
    :param: nval the index of the function positive integer
     :param: the xvalue at which the function should be evaluated"""

    norm_cst = 1 / (np.sqrt(2 ** (nval) * fact(nval) * np.sqrt(np.pi)))
    return norm_cst * np.exp((-(xval**2)) / 2) * hermite(nval, xval)


x_list = np.linspace(-4, 4, 100)  # list of x values

for n in range(0, 4):
    function_val_list = []
    # initialize list which will
    # contain values of harmonic oscillator func and reinitialize after each loop over n
    for x in x_list:
        function_val = harmoscillator(n, x)
        function_val_list.append(function_val)
    plot(x_list, function_val_list, label="n = " + str(n))
    title("Harmonic oscillator functions")
    legend()
    xlabel("position")
    ylabel("wavefunction")
show()

# ==================================================================================================
# Part b
# Performing the same calculation for n =30 from -10 t0 10

# x2_list = np.linspace(-10,10,50) #x values
# funcvals = harmoscillator(30,x2_list) #calculating the function

# plotting the function
# plot(x2_list,funcvals,label="n = 30 ")
# legend()
# show()

# ==================================================================================================
# Part c
# Determining the position expectation value in the harmonic oscillator
# Need to define a function which does the change of variable for us to perform the integral over
# finite ranges


def changevar(zval):
    """This function will perform the change of variable to have an integrand that
    can be evaluated on finite range. The sub is done such that -infty to infty becomes -1 to 1
    :param: zval the value at which the integral is to be evaluated"""

    x_param = zval / (1 - zval**2)  # the actual change of variable
    coeff = (1 + zval**2) / (1 - zval**2) ** 2  # this is dx in terms of dz
    eval_func = harmoscillator(
        5, x_param
    )  # evaluating the wavefunction for n =5 at x_param
    return x_param**2 * coeff * eval_func**2  # formula of integrand


# evaluating the integral to get the square of the expectation value using gaussquadrature
expval_sq = gaussquad(changevar, -1, 1, 100)

# taking the square root to obtain the actual expectation value
exp_val = np.sqrt(expval_sq)

# printing
print(exp_val)
