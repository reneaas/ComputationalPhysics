import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from linjetilpasning import Linjetilpasning
from matplotlib import rc
rc('text', usetex=True)

filename = "results/BucklingBeam/benchmarks/n_vs_iterations/n_vs_iterations_collected.txt"
N = []
Iterations = []
with open(filename, "r") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        numbers = line.split()
        N.append(float(numbers[0]))
        Iterations.append(float(numbers[-1]))

Iterations = [np.log10(iterations) for iterations in Iterations]
N = [np.log10(n) for n in N]

N = np.array(N)
Iterations = np.array(Iterations)
Line = Linjetilpasning(N, Iterations)
slope, dslope, I0, dI0 = Line.linjetilpasning()
print(slope)

figurename = "plot_n_vs_iterations.pdf"
path = "results/BucklingBeam/benchmarks/plots"
if not os.path.exists(path):
    os.makedirs(path)

n = np.linspace(0.5, 2.5)
I = slope*n + I0

plt.scatter(N,Iterations, c = "orange", label = "datapoints")
plt.plot(n, I,c = "steelblue" ,label="slope = " + str("%.3f" % slope) + " $\pm$ " + str("%.3f" % dslope))
plt.xlabel(r"$\log_{10} (n)$" ,fontsize = 18)
plt.ylabel(r"$\log_{10} (I)$", fontsize = 18)
plt.legend(fontsize = 14)
plt.xticks(size = 12)
plt.yticks(size = 12)
plt.savefig(figurename, dpi = 1000)

os.system("mv" + " " + figurename + " " + path)
