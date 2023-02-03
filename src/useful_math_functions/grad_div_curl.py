"""This file contains a set of function to determine the gradient, divergence and curl of
scalar, and vector fields respectively"""

import numpy as np
# We shall start implementing grad, div, curl and all that haha

def gradient(function, array,stepsize):
  """This function will compute the gradient of a user specifier function"""
  x,y,z = array
  derivx = (function(np.array([x+stepsize,y,z]))-function(np.array([x-stepsize,y,z])))/(2*stepsize)
  derivy = (function(np.array([x,y+stepsize,z]))-function(np.array([x,y-stepsize,z])))/(2*stepsize)
  derivz = (function(np.array([x,y,z+stepsize]))-function(np.array([x,y,z-stepsize])))/(2*stepsize)

  return np.array([derivx,derivy,derivz])