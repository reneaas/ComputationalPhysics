import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

n = int(sys.argv[1])
outfilename = str(sys.argv[2])
angular_frequency = str(sys.argv[3])
repulsion = str(sys.argv[4])

wavefunction = []
rho = []

with open(outfilename, "r") as outfile:
    lines = outfile.readlines()
    for line in lines:
        number = line.split()
        wavefunction.append(float(number[0]))
        rho.append(float(number[-1]))

if repulsion == "yes":
    path = "results/QM_TwoElectrons/Repulsion/plots"
    figurename = "groundstate_n_" + str(n) + "_omega_"  + angular_frequency + "_with_repulsion" + ".png"
else:
    path = "results/QM_TwoElectrons/NoRepulsion/plots"
    figurename = "groundstate_n_" + str(n) + "_omega_"  + angular_frequency + "_no_repulsion"  +  ".png"

myfontsize = 18
plt.plot(rho, wavefunction, label = "Ground state electron")
plt.xlabel(r"$\rho$", fontsize = myfontsize)
plt.ylabel(r"$|\psi ( \rho )|^2$", fontsize = myfontsize)
plt.savefig(figurename, dpi = 1000)

if not os.path.exists(path):
    os.makedirs(path)
os.system("mv" + " " + figurename + " " + path)
