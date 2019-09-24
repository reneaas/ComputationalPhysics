import numpy as np
class Linjetilpasning:
    """
    How to use this class:
    Suppose we have to arrays x and y where y(x) = A*x + B.
    An instance is then made as

    Line = Linjetilpasning(x,y)

    To compute the constanst A and B and their corresponding
    standard deviations dA and dB, write the following:

    A, B, dA, dB = Line.linjetilpasning()

    """
    def __init__(self, x, y, number_of_datasets):
        self.x = x
        self.y = y
        self.number_of_datasets = number_of_datasets
        if self.number_of_datasets != 1:


    def linjetilpasning(self):
        #linjetilpasning y = mx + c med usikkerheter dm og dc
        x = self.x
        y = self.y
        n = len(x)

        #algoritme
        D = np.dot(x,x) - (1/n)*np.sum(x)**2
        E = np.dot(x,y) - (1/n)*np.sum(x)*np.sum(y)
        F = np.dot(y,y) - (1/n)*np.sum(y)**2
        x_mean = (1/n)*np.sum(x)
        y_mean = (1/n)*np.sum(y)
        m = E/D
        self.c = y_mean - m*x_mean
        dm_squared = (1/(n-2))*(D*F-E**2)/D**2
        dc_squared = (1/(n-2))*(D/n + x_mean**2)*(D*F-E**2)/D**2
        self.dm = np.sqrt(dm_squared)
        self.dc = np.sqrt(dc_squared)
        return m, c, dm, dc

    def make_plot(self, x_datapoints, y_datapoints):
        self.X = x_datapoints
        self.Y = y_datapoints
