"""In this file we will write a function to determine the acceleration
and velocity of an object from a set of coordinates."""
import numpy as np
import matplotlib.pyplot as plt


def velocity(data):
    """
    This function will determine the velocity from position vs time data.
    The method of choice will be the central difference method for all except
    first and last point, for which we'll use forward and backward difference
    respectively (this is not great and will lead to errors in acceleration at these points)
    the user may choose to omit the first and last points
    The function will return an array of velocity values and
    an array of time values. This is done such that we can call the function
    on itself to get acceleration array later on

    param: data: which is the data set

    In this case, the data set is understood as the transpose of a numpy array
    This function does not make use of loops, instead we will generate all
    derivative values at once and store them into an array
    """

    velocityarray = np.zeros(len(data[:, 1]))
    # initializing an array of zeros the length of the y-data array

    velocityarray[1:-1] = (data[2:, 1] - data[:-2, 1]) / (data[2:, 0] - data[:-2, 0])
    # setting up the central difference method for all points except first and last
    # Essentially we are replacing all zeros by the value of the derivative

    # setting up forward difference method for first point
    velocityarray[0] = (data[1, 1] - data[0, 1]) / (data[1, 0] - data[0, 0])

    # setting up backward difference method for last point
    velocityarray[-1] = (data[-1, 1] - data[-2, 1]) / (data[-1, 0] - data[-2, 0])
    return np.transpose((data[:, 0], velocityarray))


if __name__ == "__main__":
    # loading the data file and storing it into a variable that contains position vs time
    # only works this way if the file is in the same directory as the python script
    # knowing the structure of your data file is important before you start processing it
    initdata = np.loadtxt("04b_Exercise_velocity_acceleration_data_file.dat")

    velocityvstime = velocity(initdata)
    # function call and assigning it to variable velocityvstime

    plt.scatter(velocityvstime[:, 0], velocityvstime[:, 1])
    plt.xlabel("times [s]")
    plt.ylabel("Velocity [m/s]")
    plt.show()

    # calculating acceleration using central difference
    accelerationvstime = velocity(velocity(initdata))

    plt.plot(accelerationvstime[:, 0], accelerationvstime[:, 1], "r+")
    plt.xlabel("time [s]")
    plt.ylabel("Acceleration [m/s^2]")
    plt.show()
