"""In this file, we will write codes to compute integrals over infinite ranges"""

import numpy as np
from scipy.constants import c as _c, Boltzmann as _kb, hbar, Stefan_Boltzmann as sigma

from src.useful_math_functions.integrals import gaussquad

# ================================================================================================
# Stefan-Boltzmann constant
# taken from exercise 5.12 (p181) from Newman Computational Physics

# We will evaluate the integral of the total rate at which energy is radiate
# Recall math W = (k_b^4 T^4)/(4 pi^2 c^2 h_bar^3) integral_0 to inf (x^3/(e^x -1)


def integrand(zval):
    """This function simply calculates the integrand needed for the integral to determine the
    total rate at which energy is radiated
    :param: zval which is the dependent variable of the function
    """
    return ((zval**3) / (1 - zval) ** 5) * (1 / (np.exp(zval / (1 - zval)) - 1))


def energyrate(temperature):
    """Caclculates the total energy by performing the integration of integrand function
    for a given temperature, it uses gaussian quadrature to do so
    :param: TEMPERATURE which is a constant
    """

    pre_factor = (_kb**4 * temperature**4) / (4 * np.pi**2 * _c**2 * hbar**3)
    integral = gaussquad(integrand, 0, 1, 50)
    total_energyrate = pre_factor * integral
    return total_energyrate
# ================================================================================================


if __name__ == "__main__":
    # Setting the value of the temperature
    TEMPERATURE = 293  # in Kelvin

    # Calculating the total rate
    rate = energyrate(TEMPERATURE)

    # Calculating the Stefan_Botlzmann constant and the error
    stef_boltz = rate / TEMPERATURE**4
    error = abs(stef_boltz - sigma)
    print(
        "The total rate is ",
        rate,
        "The Stefan_Boltzman constant is ",
        stef_boltz,
        "the error is ",
        error,
    )
