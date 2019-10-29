import numpy as np
import matplotlib.pyplot as plt
import sys
import os
plt.rc("text", usetex = True)

part = str(input("Which part of the project to run: [b,c,d] \n"))

if part == "c":
    T = 2.4
    path = "results/partC/"
    infilename_ordered = "MC_" + str(int(4e6)) + "_n_20_T_" + str(T) + "_ordered_.txt"
    infilename_random = "MC_" + str(int(4e6)) + "_n_20_T_" + str(T) + "_random_.txt"
    E_ordered = []
    M_ordered = []
    acceptance_ordered = []
    E_random = []
    M_random = []
    acceptance_random = []
    time = []

    with open(path + infilename_ordered, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            time.append(float(values[0]))
            E_ordered.append(float(values[1]))
            M_ordered.append(float(values[2]))
            acceptance_ordered.append(float(values[3]))

    with open(path + infilename_random, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            E_random.append(float(values[1]))
            M_random.append(float(values[2]))
            acceptance_random.append(float(values[3]))


    plt.plot(time[:], E_ordered[:], label = "ground state initiation")
    plt.plot(time[:], E_random[:], label = "random initiation")
    plt.xlabel(r"$t$ [cycles/spins]", size = 14)
    plt.ylabel(r"$\langle E \rangle $/spins", size = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 12)
    plt.figure()

    plt.plot(time[:], M_ordered[:], label = "ground state initiation")
    plt.plot(time[:], M_random[:], label = "random initiation")
    plt.xlabel("$t$ [cycles/spins]", size = 14)
    plt.ylabel(r"$\langle |M| \rangle $/spins", size = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 12)
    plt.figure()

    plt.plot(time[:], acceptance_random[:], label = "Accepted states (random initiation)")
    plt.plot(time[:], acceptance_ordered[:], label = "Accepted states (ground state initiation)")
    plt.xlabel("$t$ [cycles/spins]", size = 14)
    plt.ylabel("Accepted spins", size = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 12)
    plt.show()

if part == "d":
    energies = []
    temperature = float(input("Temperature = "))
    initial_spin_state = str(input("ordered or random: [o/r]"))
    if initial_spin_state == "o":
        initial_spin_state = "ordered"
    if initial_spin_state == "r":
        initial_spin_state = "random"

    infilename =  "boltzmann_distribution_T_" + str(temperature) + "_" + initial_spin_state + ".txt"
    path = "results/partC/"
    with open(path + infilename, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            energies.append(float(values[0]))

    plt.hist(energies, 800 + 1, density = True)
    plt.show()
