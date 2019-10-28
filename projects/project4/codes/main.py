import numpy as np
import os



print("compiling")
os.system("c++ -O3 -c main.cpp")
os.system("c++ -O3 -o main.exe main.cpp")

number_of_temperatures = int(input("Number of temperatures: "))
number_of_MC_runs = int(input("Number of Monte Carlo runs: "))



if number_of_temperatures == 1 and number_of_MC_runs == 1:
    T = float(input("Temperature: "))
    MC_cycles = int(input("Number of Monte Carlo cycles: "))
    n_spins = int(input("Number of spins in each direction = "))
    outfilename = "random.txt"
    initialize = "ordered"
    print("executing")

    command_line_args = str(number_of_temperatures) + " " + outfilename + " " + str(n_spins)\
                        + " " + str(MC_cycles) + " " + str(number_of_MC_runs) + " " + initialize + " " + str(T)

    os.system("./main.exe" + " " + command_line_args)

if number_of_temperatures == 1 and number_of_MC_runs > 1:
    T = float(input("Temperature: "))
    MC_initial= int(input("Initial Monte Carlo cycle number: "))
    MC_final = int(input("Final Monte Carlo cycle number: "))
    n_spins = int(input("Number of spins in each direction = "))
    MC_stepsize = float(MC_final-MC_initial)/float(number_of_MC_runs)
    outfilename = "thermodynamic_quantities_T_" + str(T) + ".txt"
    initialize = "ordered"
    MC_cycles = 1
    outfilename2 = "Relative_error_n_2.txt"
    print("executing")

    command_line_args = str(number_of_temperatures) + " " + outfilename + " " + str(n_spins) \
                        + " " + str(MC_cycles) + " " + str(number_of_MC_runs) + " " + initialize + " " + str(T) \
                        + " " + str(MC_initial) + " " + str(MC_stepsize)

    os.system("./main.exe" + " " + command_line_args)
    files = outfilename + " " + outfilename2
    path = "results/2x2"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + files + " " + path)
    print("Finito.")

if number_of_temperatures > 1:
    outfilename = "thermodynamic_quantities.txt"
    T_initial = float(input("T_initial? "))
    T_final = float(input("T_final? "))
    step_size = (T_final - T_initial)/number_of_temperatures
    MC_cycles = int(input("Number of Monte Carlo cycles: "))
    n_spins = int(input("Number of spins in each direction = "))
    initialize = "ordered"
    number_of_MC_runs = "hei"
    command_line_args = str(number_of_temperatures) + " " + outfilename + " " + str(n_spins)\
                        + " " + str(MC_cycles) + " " + number_of_MC_runs + " " + initialize + " " + str(T_initial) + " " + str(T_final) + " " + str(step_size)

    os.system("./main.exe" + " " + command_line_args)
