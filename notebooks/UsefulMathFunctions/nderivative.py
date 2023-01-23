"""This file contains a function to approximate the nth-derivative of functions."""
from math import factorial as fact


def nderivative(function, order, argument, stepsize):
    """This function outputs the nth derivative
    of a user-specified function.

    It takes as input the following
    :param function: a user-defined function
    :param order: the order of the derivative of interests
    :param argument: the value at which the derivative is to be evaluated
    :param stepsize: the size of the intervals: warning: to small a value leads to large issues
    Author : Kedy Edme
    Date : Jan 20th 2023
    """
    # declaring the derivative variable and setting it to zero
    derivative = 0

    # setting up a loop to perform the sum for the nth derivative

    for k in range(order + 1):
        derivative = derivative + (
            (-1) ** (order + k)
            * (fact(order) / (fact(k) * fact(order - k)))
            * function(argument + k * stepsize)
        ) / (stepsize**order)

    return derivative
