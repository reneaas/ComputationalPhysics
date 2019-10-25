import numpy as np
import os



print("compiling")
os.system("c++ -O3 -c main.cpp")
os.system("c++ -O3 -o main.exe main.cpp")
print("executing")
os.system("./main.exe lol 20 100 1")


number_of_temperatures = int(input("Number of temperatures: "))

if number_of_temperatures == 1:
    T = float(input("T: "))
else:

    T_initial = float(input("T_initial? "))
    T_final = float(input("T_final? "))
    step_size = (T_final - T_initial)/number_of_temperatures


"""
Rekkefølge på shit c++ programmet tar:
antall temperaturer - outfilename - antall dimensjoner på matrise - antall monte carlo samples - T_initial - T_final - step_size
"""
