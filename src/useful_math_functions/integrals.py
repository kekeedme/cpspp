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


def adaptivetz(function, avalue, bvalue, nslices, target_error):
    """This function will perform the integration using adaptive trapezoid method
    it uses
    :param: function the function we are integrating
    :param: avalue and bvalue the upper and lower bound of the interval of interest
    :param: nslices intial number of slices N
    :param: target_error the target error we wish to achieve
    It returns the value of the integral AND the calculated error
    N.B: This function should be imported with the function_integrationtz() function"""

    # Recall math: (0.5)*I_i_less1 + h * sum_odd(f(a+kh))

    # calculating the first integral using the trapezoid rule (so the I_i_less1 term)
    integral_1 = function_integrationtz(function, avalue, bvalue, nslices)
    integral_new = None  # will be used and updated later in the while loop
    actual_error = (1 / 3) * (integral_1 - 0)  # defining a VERY approximate
    # initial value of the error to be used in the while loop since we have not yet the second I_i

    while abs(actual_error) > target_error:
        nslices *= 2  # double the number of slices if we have not reached the target
        width = (bvalue - avalue) / nslices  # redefine the width for each new slice
        integral_n = function_integrationtz(
            function, avalue, bvalue, nslices
        )  # recalculate the ith-1 integral
        # generating the sum over odd terms
        sum_adapt = 0
        for k in range(1, nslices, 2):
            sum_adapt += function(avalue + k * width)

        integral_new = (0.5 * integral_n) + (
            width * sum_adapt
        )  # calculate the ith integral using Ii = (1/2)I_i-1 + hi * sum_odd (f(a+khi))
        actual_error = (1 / 3) * (
            integral_new - integral_n
        )  # recalculate the error from the previous integrals
    return integral_new, actual_error


def Rombergtz(function, avalue, bvalue, nslices, target_error):
    """This function will perform the integration using Romberg trapezoid method
    it uses
    :param: function the function we are integrating
    :param: avalue and bvalue the upper and lower bound of the interval of interest
    :param: nslices intial number of slices N
    :param: target_error the target error we wish to achieve
    It returns the value of the integral AND the calculated error
    N.B: This function should be imported with the function_integrationtz() function"""

    # calculate the integral for the first the very first integral and set it equal to integral_1
    integral_1_1 = function_integrationtz(function, avalue, bvalue, nslices)

    m_index = 1  # index counting the number of terms in Romberg
    actual_error = 1  # defining an approximate error
    inte_list_i_less = (
        []
    )  # creating a list to store all necessary Romberg calcs and use them for further calculations
    inte_list_i_less.append(integral_1_1)  # Add first value in list: R1,1
    while abs(actual_error) > target_error:
        nslices *= 2  # update the number of slices by multiplying by two
        m_index += 1  # increase the Romberg index
        inte_list_i = (
            []
        )  # make a new list that will contain the latest values of the Romberg terms
        integral_n = function_integrationtz(
            function, avalue, bvalue, nslices
        )  # Calculate next integral with twice as many slices as before
        inte_list_i.append(integral_n)  # add to that new list
        # loop that will carry out the calc to get higher Romberg terms from previous values of lists, since we already did m=1 we need to start at m=2
        for m in range(2, m_index + 1):
            # Calculating the R_i_m+1 term
            r_i_m_plus_one = inte_list_i[-1] + (
                (1 / (4 ** (m) - 1)) * (inte_list_i[-1] - inte_list_i_less[m - 2])
            )  # calculate Romberb Ri,m_plus1
            inte_list_i.append(
                r_i_m_plus_one
            )  # append this value since it will be needed for subsequent Romberg terms
            actual_error = (1 / (4**m - 1)) * (
                inte_list_i[-1] - inte_list_i_less[m - 2]
            )  # update the error on the penultimate value
        inte_list_i_less = (
            inte_list_i  # update the first list to become the latest list
        )

    return inte_list_i_less[-1], actual_error


def gaussquad(function, avalue, bvalue, npoints):
    """This function performs integration using the Gaussian quadrature rule
    :param: function the function we which to integrate
    :param: avalue and bvalue the upper and lower bound of the interval of interest
    :param: npoints which is the number of points required remember polynomial
    need to be of degree 2N-1"""
    from gaussxw import gaussxw

    # needed from the gaussxw that calculates the weight and points from Gauss-Legendre
    # for nth Legendre function. They are ARRAYS
    x_points, weights = gaussxw(npoints)
    # Remapping the array of points to suit our interval and therefore rescaling weights as well
    mapped_x = (0.5) * (bvalue - avalue) * x_points + 0.5 * (bvalue + avalue)
    mapped_weights = 0.5 * (bvalue - avalue) * weights

    # Calculate the integral
    integral = 0.0
    for k in range(npoints):
        integral += mapped_weights[k] * function(mapped_x[k])
    return integral
