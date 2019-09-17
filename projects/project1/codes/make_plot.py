import numpy as np
import matplotlib.pyplot as plt
import os
import sys

#filename  = str(input("Type filename with extension: "))
filename = str(sys.argv[1])
which_algorithm = str(sys.argv[2])
v = []
with open(filename, "r") as infile:
    lines = infile.readlines()
    for line in lines:
        numbers = line.split()
        for number in numbers:
            v.append(float(number))

v = np.array(v)

def u(x):
    """
    Analytical solution to the differential equation
    """
    return 1 - (1-np.exp(-10))*x - np.exp(-10*x)

n = len(v)
h = 1/(float(n) + 1)
x = [float((i+1)*h) for i in range(n)]
x = np.array(x)

plt.plot(x,v, "-b",label="Numerical solution with n = " + str(n))
plt.plot(x,u(x), "--r",  label="Analytical solution")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()

figurename = filename.strip(".txt") + ".png"
plt.savefig(figurename)
if which_algorithm == "general_algorithm":
    directory = "results/general_algorithm/plots/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.system("mv " + figurename + " " + directory)
if which_algorithm == "special_algorithm":
    directory = "results/special_algorithm/plots/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.system("mv " + figurename + " " + directory)
if which_algorithm == "LU":
    directory = "results/LU/plots/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.system("mv " + figurename + " " + directory)
plt.close()
