import numpy as np
import matplotlib.pyplot as plt

class StraightLine:
    """
    How to use this class:

    ----------------------Step 1: Initiate an instance--------------------------

    Lines = StraightLine(x, y, number_of_datasets)

    x: (n x m) - matrix with n datasets of length m
    y: (n x m) - matrix with n datasets of length m
    number_of_datasets: the number of datasets n.

    -------------------Step 2: Find the straight lines.------------------------

    To find the straight line approximations and standard deviations for all datasets,
    write

    Lines.straightline()


    ----------------Step 3: Create plots of the straight lines------------------
    This plot will include the plot of the straight lines along with its datapoints
    and the corresponding uncertainties represented as error bars. Write

    Lines.make_plot(labeltexts, figurename)

    labeltexts: a list containing all the label for each straight line.
    figurename: filename to save plot to. Recommended use: name.pdf

    """


    def __init__(self, x, y, number_of_datasets):
        self.X = x
        self.Y = y
        self.number_of_datasets = number_of_datasets

    def straightline(self):
        #linjetilpasning y = mx + c med usikkerheter dm og dc
        self.M = []
        self.C = []
        self.dM = []
        self.dC = []
        n = np.shape(self.X)[-1]
        for i in range(self.number_of_datasets):
            x = self.X[i]
            y = self.Y[i]
            #algoritme
            D = np.dot(x,x) - (1/n)*np.sum(x)**2
            E = np.dot(x,y) - (1/n)*np.sum(x)*np.sum(y)
            F = np.dot(y,y) - (1/n)*np.sum(y)**2
            x_mean = (1/n)*np.sum(x)
            y_mean = (1/n)*np.sum(y)
            m = E/D
            c = y_mean - m*x_mean
            dm_squared = (1/(n-2))*(D*F-E**2)/D**2
            dc_squared = (1/(n-2))*(D/n + x_mean**2)*(D*F-E**2)/D**2
            dm = np.sqrt(dm_squared)
            dc = np.sqrt(dc_squared)
            self.M.append(m)
            self.C.append(c)
            self.dM.append(dm)
            self.dC.append(dc)
        self.M = np.array(self.M)
        self.C = np.array(self.C)
        self.dM = np.array(self.dM)
        self.dC = np.array(self.dC)

    def make_plot(self, labeltexts, figurename):
        for i in range(self.number_of_datasets):
            labeltext = labeltexts[i]
            y_error = np.sqrt( self.dM[i]**2 + self.dC[i]**2 )
            x = np.linspace(min(self.X[i]), max(self.X[i]), 101)
            y = self.M[i]*x + self.C[i]
            plt.plot(x,y, label = labeltext)
            plt.errorbar(self.X[i], self.Y[i], y_error, capsize = 5, fmt=".")
        plt.xlabel("x", fontsize = 18)
        plt.ylabel("y", fontsize = 18)
        plt.legend(fontsize = 14)
        plt.xticks(size = 12)
        plt.yticks(size = 12)
        plt.savefig(figurename, dpi = 1000)
