import numpy as np
import os
import sys

print("__________________________________________")
print("Which part would you like to run?")
print("1-dimensional schemes          --> Type 1")
print("2-dimensional scheme           --> Type 2")
print("Stability analysis for 1D-case --> Type 3")
print("__________________________________________")

d = int(input())
methods = ["explicit", "implicit", "CN"]

if d != 22:
    print("Compiling code...")
    os.system("c++ -O3 -Wall -c main.cpp functions.cpp")
    os.system("c++ -O3 -Wall -o main.exe main.o functions.o")
    print("Compilation finished, executing program...")
else:
    print("Compiling code...")
    os.system("mpicxx -O2 -Wall -c main_mpi.cpp")
    os.system("mpicxx -O2 -o main_mpi.exe main_mpi.o")
    print("Compilation finished, executing program...")



if d == 1:   #1D case
    path = "results/1D"
    if not os.path.exists(path):
        os.makedirs(path)

    dx = [0.1, 0.01]
    total_time = [0.02, 0.085]
    r = 0.5
    for method in methods:
        print("Running method: " + method)
        for i in dx:
            if i  == 0.1:
                for t in range(len(total_time)):
                    outfilename = str(method) + "_dx_" + str(i) + "_time_" + str(t+1) + ".txt"

                    os.system("./main.exe" + " " + str(d) + " " + str(i) + " " + method + " " + outfilename + " " + str(total_time[t]) + " " + str(r))
                    os.system("mv" + " " + outfilename +" "+ path)

            if i == 0.01:
                for t in range(len(total_time)):
                    outfilename = str(method) + "_dx_" + str(i) + "_time_" + str(t+1) + ".txt"

                    os.system("./main.exe" + " " + str(d) + " " + str(i) + " " + method + " " + outfilename + " " + str(total_time[t]) + " " + str(r))
                    os.system("mv" + " " + outfilename + " " + path)


if d == 2:   #2D case
    path = "results/2D/"
    if not os.path.exists(path):
        os.makedirs(path)
    h = 0.01
    outfilename = "2D_Results.txt"
    total_time = 0.1
    r = 0.25
    os.system("./main.exe" + " " + str(d) + " " + str(h) + " " + outfilename + " " + str(total_time) + " " + str(r))
    os.system("mv" + " " + outfilename +" "+ path)

if d == 22:   #MPI dritt som ikke funker
    os.system("mpirun -np 2 ./main_mpi.exe")

if d == 3:    #Stability analysis
    path = "results/1D/Stability"
    if not os.path.exists(path):
        os.makedirs(path)

    d = "1"
    dx = 0.01
    total_time = 0.02
    r = [0.5, 0.505]
    for R in r:
        for method in methods:
            outfilename = str(method) + "_r_" + str(R) + "_tot_time_" + str(total_time) + ".txt"
            os.system("./main.exe" + " " + d + " " + str(dx) + " " + method + " " + outfilename + " " + str(total_time) + " " + str(R))
            os.system("mv" + " " + outfilename + " " + path)
