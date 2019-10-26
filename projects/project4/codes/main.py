import numpy as np
import os



print("compiling")
os.system("c++ -O3 -c main.cpp")
os.system("c++ -O3 -o main.exe main.cpp")

number_of_temperatures = int(input("Number of temperatures: "))

if number_of_temperatures == 1:
    T = float(input("T: "))
    MC_cycles = int(input("Number of Monte Carlo cycles: "))
    n_spins = int(input("Number of spins in each direction = "))
    outfilename = "random.txt"


    print("executing")
    command_line_args = str(number_of_temperatures) + " " + outfilename + " " + str(n_spins)\
                        + " " + str(MC_cycles) + " " + str(T)
    os.system("./main.exe" + " " + command_line_args)
    #antall temp - outfilename - grid - MC_cycles - temp

else:

    T_initial = float(input("T_initial? "))
    T_final = float(input("T_final? "))
    step_size = (T_final - T_initial)/number_of_temperatures
