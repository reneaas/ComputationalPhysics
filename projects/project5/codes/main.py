import numpy as np
import os
import sys


print("Compiling code...")
os.system("c++ -O3 -Wall -c main.cpp functions.cpp")
os.system("c++ -O3 -Wall -o main.exe main.o functions.o")
print("Compilation finished, executing program...")


methods = ["explicit", "implicit", "CN"]

path = "results/1D"
if not os.path.exists(path):
    os.makedirs(path)

d = 1
<<<<<<< HEAD
dx = 0.01
method = "implicit"
outfilename = "test_" + str(method) + ".txt"
=======
dx = [0.1, 0.01]

for method in methods:
    print("Running method: " + method)
    for i in dx:
        if i  == 0.1:
            time_index = [4,17]
            for t in range(len(time_index)):
                outfilename = str(method) + "_dx_" + str(i) + "_time_" + str(t+1) + ".txt"

                os.system("./main.exe" + " " + str(d) + " " + str(i) + " " + method + " " + outfilename + " " + str(time_index[t]))
                os.system("mv" + " " + outfilename +" "+ path)

        if i == 0.01:
            time_index = [10,1800]
            for t in range(len(time_index)):
                outfilename = str(method) + "_dx_" + str(i) + "_time_" + str(t+1) + ".txt"
>>>>>>> 827d2194f8dca15f91b7ea88f7bd27b3d081a506

                os.system("./main.exe" + " " + str(d) + " " + str(i) + " " + method + " " + outfilename + " " + str(time_index[t]))
                os.system("mv" + " " + outfilename + " " + path)
