import numpy as np
import os
import sys

omega = [0.01, 0.05, 0.25]
eigenvalues = []
for angular_frequency in omega:
    path = "results/QM_TwoElectrons/Repulsion/computed_eigenvalues/" \
            + "computed_eigenvalues_QM_TwoElectrons_n_100_omega_"+ \
            str(angular_frequency) + "_repulsion_yes.txt"


    with open(path, "r") as outfile:
        line = outfile.readline()
        eigenvalues.append(float(line)/2.0)

Omega = [1/angular_frequency for angular_frequency in omega]
filename = "angular_frequency_vs_gs_energy.txt"
with open(filename, "w") as infile:
    infile.write("1/omega" + " " + "energy" + "\n")
    for angular_frequency, energy in zip(Omega, eigenvalues):
        infile.write( str("%.0f" % angular_frequency) + " " + str("%.4f" % energy) + "\n" )

path = "results/QM_TwoElectrons/Repulsion"
if not os.path.exists(path):
    os.makedirs(path)

os.system("mv" + " " + filename + " " + path)

NumberOfGridpoints = [10, 20, 50, 100, 150]
outfilename = "n_vs_iterations_collected.txt"
iterations = []
for n in NumberOfGridpoints:
    filename_NumberOfIterations = "results/BucklingBeam/benchmarks/n_vs_iterations/n_vs_iterations_n_equals_" + str(n) + ".txt"
    with open(filename_NumberOfIterations, "r") as infile:
        infile.readline()
        line = infile.readline()
        numbers = line.split()
        iterations.append(float(numbers[-1]))
    #os.system("rm" + " " + filename_NumberOfIterations)


with open(outfilename, "w") as outfile:
    outfile.write("n" + " " + "iterations" + "\n")
    for n, iter in zip(NumberOfGridpoints, iterations):
        outfile.write(str(n) + " " + str(iter) + "\n")

path_n_vs_iterations = "results/BucklingBeam/benchmarks/n_vs_iterations"
os.system("mv" + " " + outfilename + " " + path_n_vs_iterations)
