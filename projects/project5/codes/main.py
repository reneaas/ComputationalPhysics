import numpy as np
import os
import sys



print("Compiling code...")
os.system("c++ -O3 -Wall -c main.cpp functions.cpp")
os.system("c++ -O3 -Wall -o main.exe main.o functions.o")
print("Compilation finished, executing program...")


methods = ["explicit", "implicit", "CN"]



d = int(input("d?"))


if d == 1:
    path = "results/1D"
    if not os.path.exists(path):
        os.makedirs(path)
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

                    os.system("./main.exe" + " " + str(d) + " " + str(i) + " " + method + " " + outfilename + " " + str(time_index[t]))
                    os.system("mv" + " " + outfilename + " " + path)


if d == 2:
    path = "results/2D/"
    if not os.path.exists(path):
        os.makedirs(path)
    h = 0.01
    outfilename = "2D_Results.txt"
    time_index = 300
    os.system("./main.exe" + " " + str(d) + " " + str(h) + " " + outfilename + " " + str(time_index))
    os.system("mv" + " " + outfilename +" "+ path)
