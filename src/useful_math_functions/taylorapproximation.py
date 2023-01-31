"""This file contains a function for taylor approximations."""
from math import factorial as fact
from src.useful_math_functions.nderivative import nderivative


def taylorapprox(function, argument, a_value, maxorder, stepsize):
    """This function outputs the taylor approximation
    to nth-order of a user-specified function.

    It takes:
    :param function: the user-specified function
    :param argument: the argument at which we evaluate the approximation
              usually x-values for 1-D functions, can be a list.
    :param a_value: the center around which the user seeks the approximation
    :param maxorder: the maximum order of the approximation
    :param stepsize: the stepsize used to evaluate the derivatives
    Author: Kedy Edme
    Date: Jan 20th 2023
    """
    approximation = 0

    for order in range(maxorder + 1):
        approximation = approximation + (
            nderivative(function, order, a_value, stepsize)
            * (argument - a_value) ** order
        ) / fact(order)
    return approximation
