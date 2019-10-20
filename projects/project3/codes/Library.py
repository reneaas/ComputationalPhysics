import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex = True)

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
        if self.number_of_datasets > 1:
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
        else:
            n = np.shape(self.X)[-1]
            x = self.X
            y = self.Y
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
            self.M = m
            self.C = c
            self.dM = m
            self.dC = dc

    def make_plot(self, labeltexts, xlabel, ylabel, figurename):
        if self.number_of_datasets > 1:
            errorbar_labels = ["dataset" + " " + i for i in labeltexts]
            for i in range(self.number_of_datasets):
                labeltext = labeltexts[i]
                y_error = np.sqrt( self.dM[i]**2 + self.dC[i]**2 )
                x = np.linspace(min(self.X[i]), max(self.X[i]), 101)
                y = self.M[i]*x + self.C[i]
                plt.plot(x,y, label = labeltext)
                #plt.errorbar(self.X[i], self.Y[i], y_error, capsize = 5, fmt=".", label = errorbar_labels[i])
            plt.xlabel(xlabel, fontsize = 18)
            plt.ylabel(ylabel, fontsize = 18)
            plt.legend(fontsize = 14)
            plt.xticks(size = 12)
            plt.yticks(size = 12)
            plt.savefig(figurename, dpi = 1000)
        else:
            labeltext = labeltexts
            y_error = np.sqrt( self.dM**2 + self.dC**2 )
            x = np.linspace(min(self.X), max(self.X), 101)
            y = self.M*x + self.C
            plt.plot(x,y, label = labeltext)
            plt.errorbar(self.X, self.Y, yerr = y_error, capsize = 5, fmt=".")
            plt.xlabel(xlabel, fontsize = 18)
            plt.ylabel(ylabel, fontsize = 18)
            plt.legend(fontsize = 14)
            plt.xticks(size = 12)
            plt.yticks(size = 12)
            plt.savefig(figurename, dpi = 1000)
        plt.close()


class PlottingTool:
    """
    How to use this class

    Step 1: Initiate an instance:

    Plotmaker = PlottingTool(x, y, number_of_datasets)

    x: (n x m)- matrix containing data for the x-axis, containing n datasets of length m.
    y: (n x m)-matrix containing data for the y-axis , contaning n datasets of length m.
    number_of_datasets: the number of datasets n.

    Step 2:
    Create plot simply by calling

    Plotmaker.plot(labeltexts, xlabel, ylabel, figurename, type)

    labeltexts: labels pertaining to each dataset
    xlabel: label along the x-axis
    ylabel: label along the y-axis
    figurename: filename to save the figure with
    type: plot or scatter. plot gives a regular plot, scatter gives crosses X for each point (x,y).

    """

    def __init__(self, x, y, number_of_datasets):
        self.X = np.array(x)
        self.Y = np.array(y)
        self.number_of_datasets = number_of_datasets

    def plot(self, labeltexts, xlabel, ylabel, figurename, type):
        if self.number_of_datasets == 1:
            if type == "plot":
                plt.plot(self.X, self.Y, label = labeltexts)
                plt.xlabel(xlabel, fontsize = 14)
                plt.ylabel(ylabel, fontsize = 14)
                plt.legend(fontsize = 14)
                plt.xticks(size = 14)
                plt.yticks(size = 14)
                plt.savefig(figurename, dpi = 1000)
                plt.close()

            if type == "scatter":
                plt.scatter(self.X, self.Y, marker = "x",  label = labeltexts)
                plt.xlabel(xlabel, fontsize = 18)
                plt.ylabel(ylabel, fontsize = 18)
                plt.legend(fontsize = 14)
                plt.xticks(size = 14)
                plt.yticks(size = 14)
                plt.savefig(figurename, dpi = 1000)
                plt.close()

        if self.number_of_datasets > 1:
            if type == "plot":
                for i in range(self.number_of_datasets):
                    plt.plot(self.X[i], self.Y[i], label = labeltexts[i])
                plt.xlabel(xlabel, fontsize = 18)
                plt.ylabel(ylabel, fontsize = 18)
                plt.legend(fontsize = 14)
                plt.xticks(size = 12)
                plt.yticks(size = 12)
                plt.savefig(figurename, dpi = 1000)
                plt.close()

            if type == "scatter":
                for i in range(self.number_of_datasets):
                    plt.scatter(self.X[i], self.Y[i], marker = "x" ,label = labeltexts[i])
                plt.xlabel(xlabel, fontsize = 18)
                plt.ylabel(ylabel, fontsize = 18)
                plt.legend(fontsize = 14)
                plt.xticks(size = 12)
                plt.yticks(size = 12)
                plt.savefig(figurename, dpi = 1000)
                plt.close()
