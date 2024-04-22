"""This file will also contain instances where we apply our integration function to perform
relevant physics calculations"""
import numpy as np
from src.useful_math_functions.integrals import gaussquad
from pylab import plot,show, xlabel,ylabel,xlim,ylim,title
#================================================================================================
#PERIOD of an anharmonic oscillator
#taken from exercise 5.10 (p173) from Newman Computational Physics

#We will compute the period of an anharmonic oscillator for a potential of V(x) = x^4
#a t= 0 x = a (max amplitude) and when x = 0, T= (1/4)*T the first time
def invapotential(xval):
    """This function calculate the integrand
    math:1/sqrt(V(a)-V(x))
    :param: xval only the values of x integration"""
    np.sqrt(8*MASS)* 1/(np.sqrt((amplitude**4)-(xval**4)))
#def period(amp):
    """This function calculates the period for the anharmonic oscillator
    for a given amplitude, using the potential defined above it takes
    :param: the mass of the particle
    :param: the max amplitude
    the integral goes from 0 to the user-specified amplitude"""
 #   integral = gaussquad(invapotential,0,amp,20)
 #   periodval = np.sqrt(8*MASS)*integral
 #   return periodval

if __name__ == "__main__":
    MASS = 1 #defining the mass constant
    # we will integrate for a = 0, to a = 2, and make a plot
    period_list = []
    amplitude_list = np.linspace(0.01,2,50)
    for amplitude in amplitude_list:
        period_list.append(gaussquad(invapotential,0,2,20))
    plot(amplitude_list,period_list)
    xlim(amplitude_list[0],amplitude_list[-1])
    ylim(period_list[0],period_list[-1])
    title("Period of anharmonic oscillator")
    xlabel("Amplitude")
    ylabel("Period")
    show()
