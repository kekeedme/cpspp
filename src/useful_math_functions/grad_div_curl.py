"""This file contains a set of functions to determine the gradient, divergence and curl of
scalar, and vector fields respectively"""

import numpy as np

# We shall start implementing grad, div, curl and all that :)


def gradient(function, array, stepsize):
    """This function will compute the gradient of a user specifier function
    the derivative method of choice is the central difference method
    Importantly, the function is a scalar function, even if its argument is an array
    :param function: The user-specified function
    :param array: The components of the vector-valued function up to 3-D
    :param stepsize: The stepsize with which the derivatives should be evaluated
    the function returns an array.
    """
    xval, yval, zval = array
    partial_x = (
        function(np.array([xval + stepsize, yval, zval]))
        - function(np.array([xval - stepsize, yval, zval]))
    ) / (2 * stepsize)
    partial_y = (
        function(np.array([xval, yval + stepsize, zval]))
        - function(np.array([xval, yval - stepsize, zval]))
    ) / (2 * stepsize)
    partial_z = (
        function(np.array([xval, yval, zval + stepsize]))
        - function(np.array([xval, yval, zval - stepsize]))
    ) / (2 * stepsize)

    return np.array([partial_x, partial_y, partial_z])


def divergence(function, array, stepsize):
    """This function will compute the divergence of a vector field
    Importantly, the function (vector field) is itself an array
    the function needs to be sliced approprietly
    :param function: the user-specified function
    :param array: the array at which the function is to be evaluated
    :param stepsize: The stepsize with which the derivatives should be evaluated
    the function returns a scalar.

    """
    xval, yval, zval = array
    dfunction_xdx = (
        function(np.array([xval + stepsize, yval, zval]))[0]
        - function(np.array([xval - stepsize, yval, zval]))[0]
    ) / (2 * stepsize)
    dfunction_ydy = (
        function(np.array([xval, yval + stepsize, zval]))[1]
        - function(np.array([xval, yval - stepsize, zval]))[1]
    ) / (2 * stepsize)
    dfunction_zdz = (
        function(np.array([xval, yval, zval + stepsize]))[2]
        - function(np.array([xval, yval, zval - stepsize]))[2]
    ) / (2 * stepsize)

    return dfunction_xdx + dfunction_ydy + dfunction_zdz


def curl(function, array, stepsize):
    """This function will compute the curl of a vector field
    Importantly, the function (vector field) is itself an array
    the function needs to be sliced approprietly
    :param function: the user-specified function
    :param array: the array at which the function is to be evaluated
    :param stepsize: The stepsize with which the derivatives should be evaluated
    the function returns an array.
    """
    xval, yval, zval = array
    dfunction_xdy = (
        function(np.array([xval, yval + stepsize, zval]))[0]
        - function(np.array([xval, yval - stepsize, zval]))[0]
    ) / (2 * stepsize)

    dfunction_xdz = (
        function(np.array([xval, yval, zval + stepsize]))[0]
        - function(np.array([xval, yval, zval - stepsize]))[0]
    ) / (2 * stepsize)

    dfunction_ydx = (
        function(np.array([xval + stepsize, yval, zval]))[1]
        - function(np.array([xval - stepsize, yval, zval]))[1]
    ) / (2 * stepsize)

    dfunction_ydz = (
        function(np.array([xval, yval, zval + stepsize]))[1]
        - function(np.array([xval, yval, zval - stepsize]))[1]
    ) / (2 * stepsize)

    dfunction_zdx = (
        function(np.array([xval + stepsize, yval, zval]))[2]
        - function(np.array([xval - stepsize, yval, zval]))[2]
    ) / (2 * stepsize)

    dfunction_zdy = (
        function(np.array([xval, yval + stepsize, zval]))[2]
        - function(np.array([xval, yval - stepsize, zval]))[2]
    ) / (2 * stepsize)

    return np.array(
        [
            dfunction_zdy - dfunction_ydz,
            -(dfunction_zdx - dfunction_xdz),
            dfunction_ydx - dfunction_xdy,
        ]
    )
