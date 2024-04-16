"""This file contains a series of functions to perform integral either on user-
specified functions, or on data using different rules such as trapezoid, Simpson,
adaptive-trapezoid, Romberg-trapezoid and Gauss quadrature"""

import numpy as np


def function_integrationtz(function, initp, endp, npoints):
    """This function will perform the 1-D integral of a user-define function
    using the trapezium rule, hence the tz at the end of the name
    :param: the function which we wish to integrate
    :param: npoints are in the number of points, which will be used to calculate
            the steps
    :param: initp initial point on the domain
    :param: endp point on the domain
    """
    # Recall math: I(a,b) = h[0.5*f(b)+0.5*f(a)+sum_remaining[f(a+kh)]]

    # calculating the subinterval steps = h in the mathematical formula
    steps = (endp - initp) / (npoints)

    # calculating the function value at the initial and final points
    funcvalue_init = function(initp)  # f(a)
    funcvalue_end = function(endp)  # f(b)

    # terms = all the terms in the square bracket in the math formula
    terms = 0.5 * funcvalue_init + 0.5 * funcvalue_end

    # completing the sum and returning the integral value

    for k in range(1, npoints):
        terms += function(initp + (k * steps))
    integral = steps * terms
    return integral


def errorestimator(integral1, integral2):
    """This function uses the error estimation from Trapezoid rule
    it requires the result of two integrals performed on the same function but
    where one (integral2) used twice as many number of slices as the first (integral1)
    it takes
    :param: integral1 which is the result of the integral performed on N slices
    :param: integral2 which is the result of the integral performed on 2N slices
    It returns the error estimate on the second integral (which should be the most accurate one)
    """
    return (1 / 3) * (integral2 - integral1)


def data_integrationtz(xdata, ydata):
    """This function will perform the 1-D integral of a input data
    using the trapezium rule, hence the tz at the end of the name
    :param: xdata an array of the x_values
    :param: ydata an array of the y_values
    It returns the result of the integral"""
    # Recall math: I(a,b) = h[0.5*f(b)+0.5*f(a)+sum_remaining[f(a+kh)]]

    # defining the endpoints (interval of integration)
    initp = xdata[0]
    endp = xdata[-1]

    # calculating the subinterval steps = h in the mathematical formula
    steps = (endp - initp) / (len(ydata))

    funcvalue_init = ydata[0]  # function at the end points f(a)
    funcvalue_end = ydata[-1]  # function at the end points f(b)

    # terms = all the terms in the square bracket in the math formula
    terms = 0.5 * funcvalue_init + 0.5 * funcvalue_end

    # completing the sum

    for k in range(1, len(ydata)):
        terms += ydata[k]

    integral = terms * steps
    return integral


def function_integrationsimp(function, initp, endp, npoints):
    """This function will perform the 1-D integral of a user-define function
    using Simpson's rule, hence the simp at the end of the name
    :param: npoints are in the number of points, which will be used to calculte
            the steps
    :param: initp initial point on the domain
    :param: endp point on the domain
     It returns the result of the integral"""

    # Recall math: I(a,b) = (1/3)h[f(a)+f(b)+4*sum_odd[f(a+kh)]+2*sum_even[f(a+kh)]]

    # calculating the subinterval steps = h in the mathematical formula

    steps = (endp - initp) / (npoints)

    funcvalue_init = function(initp)  # function at the end points f(a)
    funcvalue_end = function(endp)  # function at the end points f(b)

    # terms = all the terms in the square bracket in the math formula
    terms = funcvalue_init + funcvalue_end

    # calculting sum of odd terms and adding to funcvalue at end points

    for k in range(1, npoints, 2):
        terms += 4 * function(initp + (k * steps))

    # calculting sum of even terms and adding to previous sums

    for k in range(2, npoints, 2):
        terms += 2 * function(initp + (k * steps))

    integral = terms * (1 / 3) * steps
    return integral


def data_integrationsimp(xdata, ydata):
    """This function will perform the 1-D integral of a input data
    using Simpson's rule, hence the simp at the end of the name
    :param: xdata an array of the x_values
    :param: ydata an array of the y_values
      It returns the value of the integral"""

    # Recall math: I(a,b) = (1/3)h[f(a)+f(b)+4*sum_odd[f(a+kh)]+2*sum_even[f(a+kh)]]

    # calculating the subintervals
    steps = (xdata[-1] - xdata[0]) / (len(ydata))

    # setting the function values at the end equal to funcvalue_init and funcvalue_end
    funcvalue_init = ydata[0]
    funcvalue_end = ydata[-1]

    # terms = all the terms in the square bracket in the math formula
    terms = funcvalue_init + funcvalue_end

    # summing over odd terms and adding to first two terms
    for k in range(1, len(ydata), 2):
        terms += 4 * ydata[k]

    # summing over even terms and adding to the rest of the terms
    for k in range(2, len(ydata), 2):
        terms += 2 * ydata[k]

    integral = terms * (1 / 3) * steps

    return integral
