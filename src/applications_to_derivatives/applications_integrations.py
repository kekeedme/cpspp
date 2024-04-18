import numpy as np
from src.useful_math_functions.integrals import function_integrationtz
from src.useful_math_functions.integrals import adaptivetz
from src.useful_math_functions.integrals import Rombergtz
from src.useful_math_functions.integrals import gaussquad
def fonction(x):
  return (np.sin(np.sqrt(100*x)))**2

integral,error=adaptivetz(fonction,0,1,1,10**(-6))
print(integral,error)

integral_rom,error_rom=Rombergtz(fonction,0,1,1,10**(-6))
print(integral_rom,error_rom)

def function(x):
  return x**4 -2*x +1
integral_gauss=gaussquad(function,0,2,3)
print(integral_gauss)