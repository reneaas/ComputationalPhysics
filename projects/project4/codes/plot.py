import numpy as np
import matplotlib.pyplot as plt
import sys
import os

part = str(input("Which part of the project to run: [b,c,d]"))

if part == "c":
    path = "results/partC"
    infilename_ordered = "MC_" + str(int(4e7)) + "_n_20_ordered.txt"
    infilename_random = "MC_" + str(int(4e7)) + "_n_20_random.txt"
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


    plt.plot(time, E_ordered, label = "E ordered")
    plt.plot(time, E_random, label = "E random")
    plt.xlabel("t [cycles/spins]")
    plt.ylabel("E/spins")
    plt.figure()

    plt.plot(time, M_ordered, label = "M ordered")
    plt.plot(time, M_random, label = "M random")
    plt.xlabel("t [cycles/spins]")
    plt.ylabel("E/spins")

    plt.figure()
    plt.show()
