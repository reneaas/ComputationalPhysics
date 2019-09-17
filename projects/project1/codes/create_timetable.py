from os import system
import sys

N = int(sys.argv[1])
which_algorithm = str(sys.argv[2])
number_of_gridpoints = [int(10**(i)) for i in range(1,N+1)]


time_used = []
for n in number_of_gridpoints:
    if which_algorithm == "general_algorithm":
        filename_time = "timeused_" + str(n) + "_general.txt"
        with open(filename_time, "r") as infile:
            time = float(infile.readline())
            time_used.append(time)
    if which_algorithm == "special_algorithm":
        filename_time = "timeused_" + str(n) + "_special.txt"
        with open(filename_time, "r") as infile:
            time = float(infile.readline())
            time_used.append(time)
    if which_algorithm == "LU":
        filename_time = "timeused_" + str(n) + "_LU.txt"
        with open(filename_time, "r") as infile:
            time = float(infile.readline())
            time_used.append(time)


if which_algorithm == "general_algorithm":
    filename = "n_vs_time_general.txt"
    with open(filename, "w") as outfile:
        outfile.write("n" + " " + "time/s" + "\n")
        for n, time in zip(number_of_gridpoints, time_used):
            outfile.write(str("%.1E" % n) + " " + str("%f" % time) + "\n")
    system("mv" + " " + filename + " " + "~/Documents/skole/comphys/projects/project1/codes/results/general_algorithm/tables")
if which_algorithm == "special_algorithm":
    filename = "n_vs_time_special.txt"
    with open(filename, "w") as outfile:
        outfile.write("n" + " " + "time/s" + "\n")
        for n, time in zip(number_of_gridpoints, time_used):
            outfile.write(str("%.1E" % n) + " " + str("%f" % time) + "\n")
    system("mv" + " " + filename + " " + "~/Documents/skole/comphys/projects/project1/codes/results/special_algorithm/tables")
if which_algorithm == "LU":
    filename = "n_vs_time_LU.txt"
    with open(filename, "w") as outfile:
        outfile.write("n" + " " + "time/s" + "\n")
        for n, time in zip(number_of_gridpoints, time_used):
            outfile.write(str("%.1E" % n) + " " + str("%f" % time) + "\n")
    system("mv" + " " + filename + " " + "~/Documents/skole/comphys/projects/project1/codes/results/LU/tables")
