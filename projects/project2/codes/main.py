import os
import numpy as np
import matplotlib.pyplot as plt
import sys


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
    problemtype = "QM1"
    print("Solving Schrödingers eq in 3D with one electron.")
if problemtype == "qm2":
    problemtype = "QM2"
    print("Solving Schrödingers eq in 3D with two electrons.")


print("Running code for n = " + str(n))
filename = "computed_eigenvalues_" + problemtype + "_n_" + str(n) + ".txt"
os.system("./main.exe" + " " + str(n) + " " + str(max_iterations) + " " + filename + " " + problemtype)
path = "results/" + problemtype + "/computed_eigenvalues";
if not os.path.exists(path):
    os.makedirs(path)
os.system("mv" + " " +  filename + " " + path)
