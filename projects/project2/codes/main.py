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

if problemtype == "QM_OneElectron":
    rho_max = float(input("Define max. value of rho:"))
    """
    Endrer main.cpp til å ta rho_max fra terminalen

    HJELP TIL Å ENDRE KODEN med filnavnet

    """

    print("Running code for n = " + str(n))
    #filename = "computed_eigenvalues_" + problemtype + "_n_" + str(n) + ".txt"
    filename = "computed_eigenvalues_" + problemtype + "_" + str(n) + "_" + str(rho_max) + ".txt"
    os.system("./main.exe" + " " + str(n) + " " + str(max_iterations) + " " + filename + " " + problemtype + " " + str(rho_max))
    path = "results/" + problemtype + "/computed_eigenvalues";
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " +  filename + " " + path)

if problemtype == "BucklingBeam":
    print("Running code for n = " + str(n))
    filename = "computed_eigenvalues_" + problemtype + "_n_" + str(n) + ".txt"
    filename_NumberOfIterations = "n_vs_iterations_n_equals_" + str(n) + ".txt"
    filename_Time = "bb_time_" + str(n) + ".txt"
    os.system("./main.exe" + " " + str(n) + " " + str(max_iterations) + " " + filename + " " + problemtype + " " + filename_NumberOfIterations + " "  + filename_Time)
    path = "results/" + problemtype + "/computed_eigenvalues";
    path_NumberOfIterations = "results/" + problemtype + "/benchmarks" + "/n_v_iterations"
    path_time = "results/" + problemtype + "/benchmarks" + "/n_vs_time"
    if not os.path.exists(path):
        os.makedirs(path)
    if not os.path.exists(path_NumberOfIterations):
        os.makedirs(path_NumberOfIterations)
    if not os.path.exists(path_time):
        os.makedirs(path_time)
    os.system("mv" + " " +  filename + " " + path)
    os.system("mv" + " " + filename_NumberOfIterations + " " + path_NumberOfIterations)
    os.system("mv" + " " + filename_Time + " " + path_time)


if problemtype == "QM_TwoElectrons":
    print("Running code for n = " + str(n))
    repulsion = str(input("Include electron repulsion? Type yes or no: "))
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
