import numpy as np
import os


print("compiling")
os.system("c++ -O3 -Wall -c main.cpp")
os.system("c++ -O3 -Wall -o main.exe main.cpp")

part = str(input("Which part of the project would you run? [b, c, d] \n" ))

if part == "c":
    number_of_temperatures = "1"
    dimension = "20";
    MC_samples = int(input("Specify number of Monte Carlo samples: "))
    initialize_spin_matrix = "random"
    outfilename = "MC_" + str(MC_samples) + "_n_" + dimension + "_" + initialize_spin_matrix + "_.txt"
    temperature = "1"

    command_line_args = number_of_temperatures + " " + outfilename + " " + dimension \
                        + " " + str(MC_samples) + " " + initialize_spin_matrix + " " + temperature
    print("executing")
    os.system("./main.exe" + " " + command_line_args)

    path = "results/partC"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + outfilename + " " + path)
