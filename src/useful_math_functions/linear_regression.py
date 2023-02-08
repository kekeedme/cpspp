"""This file contains a simple linear regression module."""
import numpy as np
import matplotlib.pyplot as plt


def linregression(data):
    """This function performs a linear regression on a data set
    :param my_data: the numpy loaded csv file
    we will compute:
    the sum of: x values, y values, x and y values squared, and x*y values
    the function return an array containing:
    the coefficient, the intercept and correlation coefficient"""
    xvals = data[:, 0]
    yvals = data[:, 1]
    datasetlength = len(xvals)
    xyvalues = xvals * yvals
    xvalssq = xvals**2
    yvalssq = yvals**2
    sum_xvals = sum(xvals)
    sum_yvals = sum(yvals)
    sum_xyvalues = sum(xyvalues)
    sum_xvalssq = sum(xvalssq)
    sum_yvalssq = sum(yvalssq)
    # the variable assignments above are making sure we have all the quantities we need
    coefficient = ((datasetlength * sum_xyvalues) - (sum_xvals * sum_yvals)) / (
        (datasetlength * sum_xvalssq) - sum_xvals**2
    )
    # calculing the coefficient
    intercept = (sum_yvals - (coefficient * sum_xvals)) / datasetlength
    # calculating the y-intercept
    correlation = ((datasetlength * sum_xyvalues) - (sum_xvals * sum_yvals)) / (
        np.sqrt((datasetlength * sum_xvalssq) - (sum_xvals**2))
        * np.sqrt((datasetlength * sum_yvalssq) - (sum_yvals**2))
    )
    # calculting the correlation coefficient (R-squared)
    return np.array([coefficient, intercept, correlation])


if __name__ == "__main__":
    my_data = np.loadtxt("test.csv", delimiter=",", skiprows=1, usecols=np.arange(2))
    # saving the data into variable my_data
    my_regression = linregression(my_data)
    # function call and assigning output into single array my_regression
    # the next three lines will assign each values in the array to a variable
    coeff = my_regression[0]
    intercpt = my_regression[1]
    correl = my_regression[2]
    # the next two lines are global variables defined only for plotting
    xvalues = my_data[:, 0]
    yvalues = my_data[:, 1]
    # the next line will determine the y-values based on the regression, i.e the line
    line = coeff * xvalues + intercpt
    # the next few lines are for plotting, and printing the values of interests
    plt.figure(figsize=(8,5))
    plt.scatter(xvalues, yvalues)
    plt.plot(xvalues, line, "blue",label = "Regression line")
    plt.xlabel("Concentration [umol/L]")
    plt.ylabel("Absorbance [a.u.]")
    plt.legend()
    plt.show()

    print(coeff, intercpt, correl)
