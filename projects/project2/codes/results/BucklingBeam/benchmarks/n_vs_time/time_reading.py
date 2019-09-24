import matplotlib.pyplot as plt
import numpy as np
import os
from linjetilpasning import StraightLine
#import seaborn as sns

plt.rc("text", usetex = True)
n = [10, 20, 50, 100, 150]

Jacobi_time = []
Arma_time = []

for i in range(5):
    filename = "bb_time_" + str(n[i]) + ".txt"
    infile = open(filename)
    for line in infile:
        time = line.split()
        Jacobi_time.append(float(time[1]))
        Arma_time.append(float(time[3]))
    #os.system("rm" + " " + filename)

n = [np.log10(N) for N in n]

Jacobi_time = [np.log10(i) for i in Jacobi_time]
Arma_time = [np.log10(i) for i in Arma_time]


Jacobi = np.array(Jacobi_time)
Arma = np.array(Arma_time)

Y_Datasets = np.zeros((2,len(n)))
X_Datasets = np.zeros_like(Y_Datasets)
Y_Datasets[0] = Jacobi
Y_Datasets[1] = Arma
X_Datasets[0] = n
X_Datasets[1] = n

figurename = "StraightLines.pdf"
labeltexts = ["Straight line approximation; Jacobi", "Straight line approximation; Armadillo"]
number_of_datasets = np.shape(X_Datasets)[0]
Lines = StraightLine(X_Datasets, Y_Datasets, number_of_datasets)
Lines.straightline()
Lines.make_plot(labeltexts, figurename)
