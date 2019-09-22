import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


n = int(sys.argv[1])
repulsion = str(sys.argv[2])
rho_max = float(input("rho_max = "))
angular_frequencies = [0.01, 0.05, 0.1, 0.5, 1.0]

if repulsion == "yes":
    figurename = "ground_state_wavefunctions_with_repulsion_rho_max_" + str(rho_max) + ".eps"
    for angular_frequency in angular_frequencies:
        outfilename = "results/QM_TwoElectrons/Repulsion/wavefunctions/ground_state_n_"\
                    + str(n) + "_omega_" + str(angular_frequency) + "_repulsion_yes.txt"
        wavefunction = []
        rho = []
        with open(outfilename, "r") as outfile:
            lines = outfile.readlines()
            for line in lines:
                numbers = line.split()
                wavefunction.append(float(numbers[0]))
                rho.append(float(numbers[-1]))
        plt.plot(rho, wavefunction, label = r"$\omega = $" + str(angular_frequency))


    plt.ylabel(r"$|\psi (\rho )|^2$", fontsize = 18)
    plt.xlabel(r"$\rho$", fontsize = 18)
    plt.legend()
    #plt.title(r"Ground states (including repulsion) for different angular angular_frequencies")
    plt.savefig(figurename, dpi = 400)
    plt.close()

    path = "results/QM_TwoElectrons/Repulsion/plots"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + figurename + " " + path)

if repulsion == "no":
    figurename = "ground_state_wavefunctions_no_repulsion_rho_max_" + str(rho_max) + ".eps"
    for angular_frequency in angular_frequencies:
        outfilename = "results/QM_TwoElectrons/NoRepulsion/wavefunctions/ground_state_n_"\
                    + str(n) + "_omega_" + str(angular_frequency) + "_repulsion_no.txt"
        wavefunction = []
        rho = []
        with open(outfilename, "r") as outfile:
            lines = outfile.readlines()
            for line in lines:
                numbers = line.split()
                wavefunction.append(float(numbers[0]))
                rho.append(float(numbers[-1]))
        plt.plot(rho, wavefunction, label = r"$\omega = $" + str(angular_frequency))


    plt.ylabel(r"$|\psi (\rho )|^2$", fontsize = 18)
    plt.xlabel(r"$\rho$", fontsize = 18)
    plt.legend()
    #plt.title("Ground states (no repulsion) for different angular angular_frequencies")
    plt.savefig(figurename, dpi = 400)
    plt.close()

    path = "results/QM_TwoElectrons/NoRepulsion/plots"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + figurename + " " + path)
