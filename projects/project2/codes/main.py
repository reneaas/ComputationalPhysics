import os
import numpy as np
import matplotlib.pyplot as plt
import sys

compile = input("Compile anew? Type yes or no: ")
if compile == "yes":
    print("Compiling...")
    os.system("c++ -O3 -Wall -c main.cpp functions.cpp")
    print("Creating executable...")
    os.system("c++ -O3 -Wall -o main.exe main.o functions.o -larmadillo")
print("Executing...")


n = int(sys.argv[1])
max_iterations = int(sys.argv[2])
problemtype = str(sys.argv[3])

if problemtype == "bb":
    problemtype = "BucklingBeam"
    print("Solving the Buckling Beam problem...")
if problemtype == "qm1":
    problemtype = "QM_OneElectron"
    print("Solving Schrödingers eq in 3D with one electron.")
if problemtype == "qm2":
    problemtype = "QM_TwoElectrons"
    print("Solving Schrödingers eq in 3D with two electrons.")

if problemtype != "QM_TwoElectrons":
    print("Running code for n = " + str(n))
    filename = "computed_eigenvalues_" + problemtype + "_n_" + str(n) + ".txt"
    os.system("./main.exe" + " " + str(n) + " " + str(max_iterations) + " " + filename + " " + problemtype)
    path = "results/" + problemtype + "/computed_eigenvalues";
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " +  filename + " " + path)

else:
    print("Running code for n = " + str(n))
    repulsion = str(input("Include electron repulsion?, type yes or no: "))
    angular_frequency = float(input("Give the angular frequency: "))
    filename = "computed_eigenvalues_" + problemtype + "_n_" + str(n) + "_omega_" + str(angular_frequency) + "_repulsion_" + repulsion + ".txt"
    filename_wavefunction = "ground_state_n_" + str(n) + "_omega_" + str(angular_frequency) + "_repulsion_" + repulsion + ".txt"
    os.system("./main.exe" + " " + str(n) + " " + str(max_iterations) + " " + filename + " " + problemtype + " " + str(angular_frequency) + " " + filename_wavefunction + " " + repulsion)
    if repulsion == "yes":
        path = "results/" + problemtype + "/Repulsion" + "/computed_eigenvalues"
        path_wavefunction = "results/" + problemtype + "/Repulsion" + "/wavefunctions"
    else:
        path = "results/" + problemtype + "/NoRepulsion" + "/computed_eigenvalues"
        path_wavefunction = "results/" + problemtype + "/NoRepulsion" + "/wavefunctions"
    if not os.path.exists(path):
        os.makedirs(path)
    if not os.path.exists(path_wavefunction):
        os.makedirs(path_wavefunction)

    os.system("python3 make_plot.py" + " " + str(n) + " " + filename_wavefunction + " " + str(angular_frequency) + " " + repulsion)              #Makes a plot of the wavefunction
    #Moves the written files to a their respective destinations.
    os.system("mv" + " " +  filename + " " + path)
    os.system("mv" + " " +  filename_wavefunction + " " + path_wavefunction)
