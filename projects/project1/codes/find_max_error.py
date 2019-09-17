from os import system
import sys
import numpy as np
import matplotlib.pyplot as plt

#input from terminal
N = int(sys.argv[1])

which_algorithm = str(sys.argv[2])


number_of_gridpoints = [int(10**(i)) for i in range(1,N+1)]

max_error = []
if which_algorithm == "general_algorithm":
    for n in number_of_gridpoints:
        filename_errors = "errors_n_" + str(n) + "_general.txt"
        errors = []
        with open(filename_errors, "r") as infile:
            #timeused = float(infile.readline())   #Time used to run the c++ program.
            lines = infile.readlines()
            for line in lines:
                numbers = line.split()
                for number in numbers:
                    errors.append(float(number))
        max_error.append(max(errors))

if which_algorithm == "special_algorithm":
    for n in number_of_gridpoints:
        filename_errors = "errors_n_" + str(n) + "_special.txt"
        errors = []
        with open(filename_errors, "r") as infile:
            #timeused = float(infile.readline())   #Time used to run the c++ program.
            lines = infile.readlines()
            for line in lines:
                numbers = line.split()
                for number in numbers:
                    errors.append(float(number))
        max_error.append(max(errors))

if which_algorithm == "LU":
    for n in number_of_gridpoints:
        filename_errors = "errors_n_" + str(n) + "_LU.txt"
        errors = []
        with open(filename_errors, "r") as infile:
            #timeused = float(infile.readline())   #Time used to run the c++ program.
            lines = infile.readlines()
            for line in lines:
                numbers = line.split()
                for number in numbers:
                    errors.append(float(number))
        max_error.append(max(errors))


#If the algorithm used is general Thomas
if which_algorithm == "general_algorithm":
    outfilename = "max_errors_general.txt"
    with open(outfilename, "w") as outfile:
        outfile.write("n" + " " + "max-error" + "\n")
        for i,e in zip(number_of_gridpoints, max_error):
            outfile.write(str("%.1E" % i) + " " + str("%f" % e) + "\n")

    table_directory = "~/Documents/skole/comphys/projects/project1/codes/results/general_algorithm/tables"
    system("mv " + outfilename + " " + table_directory)

    h = [1/(float(i) + 1) for i in number_of_gridpoints]
    log_h = [np.log10(i) for i in h]

    title = "Max relative error using the general Thomas algorithm"
    figurename = "max_error_general.png"
    plt.plot(log_h, max_error)
    plt.xlabel("log10(h)")
    plt.ylabel("log10(max error)")
    plt.title(title)
    plt.savefig(figurename)
    plot_directory = "~/Documents/skole/comphys/projects/project1/codes/results/general_algorithm/plots"
    system("mv " + figurename + " " + plot_directory)
    plt.close()


#if the algorithm used in Special Thomas
if which_algorithm == "special_algorithm":
    outfilename = "max_errors_special.txt"
    with open(outfilename, "w") as outfile:
        outfile.write("n" + " " + "max-error" + "\n")
        for i,e in zip(number_of_gridpoints, max_error):
            outfile.write(str("%.1E" % i) + " " + str("%f" % e) + "\n")

    table_directory = "~/Documents/skole/comphys/projects/project1/codes/results/special_algorithm/tables"
    system("mv " + outfilename + " " + table_directory)

    h = [1/(float(i) + 1) for i in number_of_gridpoints]
    log_h = [np.log10(i) for i in h]

    title = "Max relative error using the specialized Thomas algorithm"
    figurename = "max_error_special.png"
    plt.plot(log_h, max_error)
    plt.xlabel("log10(h)")
    plt.ylabel("log10(max error)")
    plt.title(title)
    plt.savefig(figurename)
    plot_directory = "~/Documents/skole/comphys/projects/project1/codes/results/special_algorithm/plots"
    system("mv " + figurename + " " + plot_directory)
    plt.close()

#if the algorithm used is LU-Thomas
if which_algorithm == "LU":
    outfilename = "max_errors_LU.txt"
    with open(outfilename, "w") as outfile:
        outfile.write("n" + " " + "max-error" + "\n")
        for i,e in zip(number_of_gridpoints, max_error):
            outfile.write(str("%.1E" % i) + " " + str("%f" % e) + "\n")

    table_directory = "~/Documents/skole/comphys/projects/project1/codes/results/LU/tables"
    system("mv " + outfilename + " " + table_directory)

    h = [1/(float(i) + 1) for i in number_of_gridpoints]
    log_h = [np.log10(i) for i in h]

    figurename = "max_error_LU.png"
    title = "Max relative error using LU-decomposition with the Thomas algorithm"
    plt.plot(log_h, max_error)
    plt.xlabel("log10(h)")
    plt.ylabel("log10(max error)")
    plt.title(title)
    plt.savefig(figurename)
    plot_directory = "~/Documents/skole/comphys/projects/project1/codes/results/LU/plots"
    system("mv " + figurename + " " + plot_directory)
    plt.close()
