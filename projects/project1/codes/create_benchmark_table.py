
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
N = int(sys.argv[1])

general_directory = "results/general_algorithm/tables/"
general_table = "n_vs_time_general.txt"

special_directory = "results/special_algorithm/tables/"
special_table = special_directory + "n_vs_time_special.txt"

LU_directory = "results/LU/tables/"
LU_table = LU_directory + "n_vs_time_LU.txt"

general_time = []
special_time = []
LU_time = []
number_of_gridpoints = [int(10**i) for i in range(1,N+1)]

filepath = os.path.join(general_directory, general_table)
infile = open(filepath)
infile.readline()
lines = infile.readlines()
for line in lines:
    numbers = line.split()
    general_time.append(float(numbers[-1]))
infile.close()

with open(special_table, "r") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        numbers = line.split()
        special_time.append(float(numbers[-1]))

with open(LU_table, "r") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        numbers = line.split()
        LU_time.append(float(numbers[-1]))

outfilename = "n_vs_time_all_algorithms.txt"
with open(outfilename, "w") as outfile:
    outfile.write("n" + " " + "general" + " " + "special" + " " + "LU" + "\n")
    for n, general, special, lu in zip(number_of_gridpoints, general_time, special_time, LU_time):
        outfile.write(str("%.1E" % n) + " " + str("%f" % general) + " " +  str("%f" % special) + " " +  str("%f" % lu) + "\n")

results_directory = "~/Documents/skole/comphys/projects/project1/codes/results/"
os.system("mv" + " " + outfilename + " " + results_directory)
