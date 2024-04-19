"""This file of functions written by Mark Newman calculates the roots to the Nth
 legendre polynomial to be used in Gauss quadrature"""
######################################################################
#
# Functions to calculate integration points and weights for Gaussian
# quadrature
#
# x,w = gaussxw(N) returns integration points x and integration
#           weights w such that sum_i w[i]*f(x[i]) is the Nth-order
#           Gaussian approximation to the integral int_{-1}^1 f(x) dx
# x,w = gaussxwab(N,a,b) returns integration points and weights
#           mapped to the interval [a,b], so that sum_i w[i]*f(x[i])
#           is the Nth-order Gaussian approximation to the integral
#           int_a^b f(x) dx
#
# This code finds the zeros of the nth Legendre polynomial using
# Newton's method, starting from the approximation given in Abramowitz
# and Stegun 22.16.6.  The Legendre polynomial itself is evaluated
# using the recurrence relation given in Abramowitz and Stegun
# 22.7.10.  The function has been checked against other sources for
# values of N up to 1000.  It is compatible with version 2 and version
# 3 of Python.
#
# Written by Mark Newman <mejn@umich.edu>, June 4, 2011
# You may use, share, or modify this file freely
#
######################################################################

from numpy import ones,copy,cos,tan,pi,linspace

def gaussxw(npoints):
    """This function calculates the roots of the Nth legendre polynomial"""
    # Initial approximation to roots of the Legendre polynomial
    avalue = linspace(3,4*npoints-1,npoints)/(4*npoints+2)
    xvalue = cos(pi*avalue+1/(8*npoints*npoints*tan(avalue)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0array = ones(npoints,float)
        p1array = copy(xvalue)
        for k in range(1,npoints):
            p0array,p1array = p1array,((2*k+1)*xvalue*p1array-k*p0array)/(k+1)
        dp = (npoints+1)*(p0array-xvalue*p1array)/(1-xvalue*xvalue)
        dx = p1array/dp
        xvalue -= dx
        delta = max(abs(dx))

    # Calculate the weights
    weights = 2*(npoints+1)*(npoints+1)/(npoints*npoints*(1-xvalue*xvalue)*dp*dp)

    return xvalue,weights

def gaussxwab(npoints,aval,bval):
    """This function remaps the x points and the weights in quadrature"""
    xvals,weights = gaussxw(npoints)
    return 0.5*(bval-aval)*xvals+0.5*(bval+aval),0.5*(bval-aval)*weights
