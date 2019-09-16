import os
import numpy as np
import matplotlib.pyplot as plt
import sys


print("Compiling...")
os.system("c++ -O3 -c main.cpp functions.cpp -larmadillo")
print("Creating executable...")
os.system("c++ -O3 -o main.exe main.o functions.o -larmadillo")
print("Executing...")


N = int(sys.argv[1])
max_iterations = int(sys.argv[2])
Number_of_Gridpoints = [int(i) for i in range(3,N)]
for n in Number_of_Gridpoints:
    filename = "computed_eigenvalues_" + "n_" + str(n) + ".txt"
    os.system("./main.exe" + " " + str(n) + " " + str(max_iterations) + " " + filename)
    path = "results/computed_eigenvalues";
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " +  filename + " " + path)
