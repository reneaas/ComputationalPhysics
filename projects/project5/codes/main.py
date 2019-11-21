import numpy as np
import os


print("Compiling code...")
os.system("c++ -O3 -Wall -c main.cpp functions.cpp")
os.system("c++ -O3 -Wall -o main.exe main.o functions.o")
print("Compilation finished, executing program...")


d = 1
dx = 0.01
method = "implicit"
outfilename = "test.txt"

os.system("./main.exe" + " " + str(d) + " " + str(dx) + " " + method + " " + outfilename)
