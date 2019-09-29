import os
import numpy as np
import matplotlib.pyplot as plt
import sys

#compile = input("Compile anew? Type yes or no: ")
compile = str(sys.argv[4])
if compile == "yes":
    print("Compiling...")
    os.system("c++ -O3 -Wall -c main.cpp functions.cpp")
    print("Creating executable...")
    os.system("c++ -O3 -Wall -o main.exe main.o functions.o -larmadillo")
print("Executing...")



n = int(sys.argv[1])
max_iterations = int(sys.argv[2])
problemtype = str(sys.argv[3])

rho_maxes = [2.0,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,\
                5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,8.0,9.0,10.0,11.0,\
                12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,25.0,30.0,35.0,40.0,45.0,50.0,100.0,200.0]



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

    run_several_iterations = str(input("Run code for several infinities? Type yes or no: "))
    if run_several_iterations == "yes":
        for rho_max in rho_maxes:
            print("Running code for rho_max = " + str(rho_max))
            #filename = "computed_eigenvalues_" + problemtype + "_n_" + str(n) + ".txt"
            filename = "computed_eigenvalues_" + problemtype + "_" + str(n) + "_" + str(rho_max) + ".txt"
            os.system("./main.exe" + " " + str(n) + " " + str(max_iterations) + " " + filename + " " + problemtype + " " + str(rho_max))
            path = "results/" + problemtype + "/computed_eigenvalues";
            if not os.path.exists(path):
                os.makedirs(path)
            os.system("mv" + " " +  filename + " " + path)
        os.system("python3 error_eigen.py")

    if run_several_iterations == "no":
        rho_max = float(input("Define max. value of rho:"))
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
