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
slope, I0, dslope, dI0 = Line.linjetilpasning()
print(slope)

figurename = "plot_n_vs_iterations.pdf"
path = "results/BucklingBeam/benchmarks/plots"
if not os.path.exists(path):
    os.makedirs(path)

n = np.linspace(0.5, 2.5)
I = slope*n + I0
print("SDV I = ", dI0)

y_error = np.sqrt(dI0**2 + dslope**2)*Iterations
#y_error = np.sqrt(dI0**2 + dslope**2)

plt.plot(n, I,c = "k" ,label="$\log_{10}(I) =$"  + " " + str("%.3f" % slope) + "$n$" + " $\pm$ " + str("%.3f" % dslope))
plt.errorbar(N, Iterations, yerr = y_error, label='Datapoints', capsize = 5, fmt = ".r")
plt.xlabel(r"$\log_{10} (n)$" ,fontsize = 18)
plt.ylabel(r"$\log_{10} (I)$", fontsize = 18)
plt.legend(fontsize = 14)
plt.xticks(size = 12)
plt.yticks(size = 12)
plt.savefig(figurename, dpi = 1000)

os.system("mv" + " " + figurename + " " + path)
