import numpy as np
import os
import sys


print("Compiling code...")
os.system("c++ -O3 -Wall -c main.cpp functions.cpp")
os.system("c++ -O3 -Wall -o main.exe main.o functions.o")
print("Compilation finished, executing program...")


method = str(input("Specify method: [explicit, implicit, CN, 2D] \n"))

if method != "2D":
    d = 1
    dx = float(input("Give dx: [0.1, 0.01] \n"))
    time_index = "100"
    outfilename = str(method) + "_dx_" + str(dx) + "_time_" + time_index + ".txt"

    os.system("./main.exe" + " " + str(d) + " " + str(dx) + " " + method + " " + outfilename + " " + time_index)
