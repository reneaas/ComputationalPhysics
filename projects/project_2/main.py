import os
import numpy as np
import matplotlib.pyplot as plt
import sys

n = int(sys.argv[1])

print("Compiling...")
os.system("c++ -O3 -c main.cpp functions.cpp -larmadillo")
print("Creating executable...")
os.system("c++ -O3 -o main.exe main.o functions.o -larmadillo")
print("Executing...")
os.system("./main.exe" + " " + str(n))
