import numpy as np
import os


print("compiling")
os.system("c++ -O3 -Wall -c main.cpp")
os.system("c++ -O3 -Wall -o main.exe main.cpp")

part = str(input("Which part of the project would you run? [b, c, d] \n" ))

if part == "b":
    spin_matrix = str(input("Would you like an ordered or random initial spin matrix? [o/r] \n" ))
    if spin_matrix == "o":
        outfilename = "Expectation_values_n_2_ordered.txt"
        outfilename2 = "Relative_error_n_2_ordered.txt"
        initialize_spin_matrix = "ordered"
    elif spin_matrix =="r":
        outfilename = "Expectation_values_n_2_random.txt"
        outfilename2 = "Relative_error_n_2_random.txt"
        initialize_spin_matrix = "random"
    else:
        print("Please choose either o or r, try again")
        os._exit(1)

    number_of_temperatures = "1"
    dimension = "2"
    MC_samples = int(input("Specify number of Monte Carlo samples: "))
    temperature = "1"


    command_line_args = number_of_temperatures + " " + outfilename + " " + dimension \
                        + " " + str(MC_samples) + " " + initialize_spin_matrix + " " + temperature + " " + outfilename2
    print("executing")
    os.system("./main.exe" + " " + command_line_args)

    path = "results/2x2"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + outfilename + " " + outfilename2 + " " + path)



if part == "c":
    number_of_temperatures = "1"
    dimension = "20";
    MC_samples = int(input("Specify number of Monte Carlo samples: "))
    initialize_spin_matrix = str(input("Ordered or randomized inital spin matrix? [o/r] \n"))
    temperature = float(input("Temperature? "))

    if initialize_spin_matrix == "o":
        initialize = "ordered"

    if initialize_spin_matrix == "r":
        initialize = "random"

    outfilename = "MC_" + str(MC_samples) + "_n_" + dimension + "_T_" + str(temperature) + "_" + initialize + "_.txt"
    outfilename2 = "boltzmann_distribution_T_" + str(temperature) + "_" + initialize + ".txt"

    command_line_args = number_of_temperatures + " " + outfilename + " " + dimension \
                    + " " + str(MC_samples) + " " + initialize + " " + str(temperature) + " " + outfilename2
    print("executing")
    os.system("./main.exe" + " " + command_line_args)

    print("moving files")
    path = "results/partC"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + outfilename + " " + outfilename2 + " " + path)
    print("Finito!!!!")
