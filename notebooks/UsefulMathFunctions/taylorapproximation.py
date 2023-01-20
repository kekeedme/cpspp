from math import factorial as fact


def nDerivative(function, order, argument, stepsize):
    """This function outputs the nth derivative
    of a user-specified function
    It takes as input the following:
    function: a user-defined function
    order: the order of the derivative of interests
    argument: the value at which the derivative is to be evaluated
    stepsize: the size of the intervals, the smaller, the better
    Author : Kedy Edme
    Date : Jan 20th 2023"""

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


def taylorapprox(function, argument, a_value, maxorder, stepsize):
    """This function outputs the taylor approximation
    to nth-order of a user-specified function. It takes:
    function: the user-specified function
    argument: the argument at which we evaluate the approximation
              usually x-values for 1-D functions, can be a list.
    a_value: the center around which the user seeks the approximation
    maxorder: the maximum order of the approximation
    stepsize: the stepsize used to evaluate the derivatives
    Author: Kedy Edme
    Date: Jan 2oth 2023"""

    approximation = 0

    for order in range(maxorder + 1):
        approximation = approximation + (
            nDerivative(function, order, a_value, stepsize)
            * (argument - a_value) ** order
        ) / fact(order)
    return approximation