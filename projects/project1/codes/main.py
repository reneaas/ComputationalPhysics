from os import system
import sys
N = int(sys.argv[1])
algorithm = str(sys.argv[2])
if algorithm == "general":
    which_algorithm = "general_algorithm"
if algorithm == "special":
    which_algorithm = "special_algorithm"
if algorithm == "lu":
    which_algorithm = "LU"

#Compile and execute program:
print("Compiling code...")
system("c++ -O3 -c project1.cpp functions.cpp")
system("c++ -O3 -o project1.exe project1.o functions.o")
print("Compilation finished, executing program...")


#run program for all n in number_of_gridpoints.
number_of_gridpoints = [int(10**(i)) for i in range(1,N+1)]
for n in number_of_gridpoints:
    if which_algorithm == "general_algorithm":
        filename_solution = "solution_part_b_n_" + str(n) + "_general.txt"
        filename_errors = "errors_n_" + str(n) + "_general.txt"
        filename_time = "timeused_" + str(n) + "_general.txt"
        print("Computing for n = " + str(n))
        system("./project1.exe" + " " + str(n) + " " + filename_solution + " " + filename_errors + " " + filename_time + " " + which_algorithm)
    if which_algorithm == "special_algorithm":
        filename_solution = "solution_part_b_n_" + str(n) + "_special.txt"
        filename_errors = "errors_n_" + str(n) + "_special.txt"
        filename_time = "timeused_" + str(n) + "_special.txt"
        print("Computing for n = " + str(n))
        system("./project1.exe" + " " + str(n) + " " + filename_solution + " " + filename_errors + " " + filename_time + " " + which_algorithm)
    if which_algorithm == "LU":
        filename_solution = "solution_part_b_n_" + str(n) + "_LU.txt"
        filename_errors = "errors_n_" + str(n) + "_LU.txt"
        filename_time = "timeused_" + str(n) + "_LU.txt"
        print("Computing for n = " + str(n))
        system("./project1.exe" + " " + str(n) + " " + filename_solution + " " + filename_errors + " " + filename_time + " " + which_algorithm)

print("Computations are done, making plots...")

#Create plots and move them to the folder ~/Documents/skole/comphys/projects/project1/codes/plots/plots_partb
for n in number_of_gridpoints:
    if which_algorithm == "general_algorithm":
        filename_solution = "solution_part_b_n_" + str(n) + "_general.txt"
        system("python3" + " " + "make_plot.py" +  " " + filename_solution + " " + which_algorithm)
        print("Plotting for n = " + str(n))
    if which_algorithm == "special_algorithm":
        filename_solution = "solution_part_b_n_" + str(n) + "_special.txt"
        system("python3" + " " + "make_plot.py" +  " " + filename_solution + " " + which_algorithm)
        print("Plotting for n = " + str(n))
    if which_algorithm == "LU":
        filename_solution = "solution_part_b_n_" + str(n) + "_LU.txt"
        system("python3" + " " + "make_plot.py" +  " " + filename_solution + " " + which_algorithm)
        print("Plotting for n = " + str(n))


print("Plots are finished, creating a plot of the maximum relative error and writes the data to a file...")


system("python3" + " " + "find_max_error.py" + " " + str(N) + " " + which_algorithm)                #calls program to find the max error, plot it and move the results to the correct directory.
system("python3" + " " + "create_timetable.py" + " " + str(N) + " " + which_algorithm)


print("Finished, removing unecessary .txt files...")
#Remove txt-files to clear up space .
for n in number_of_gridpoints:
    if which_algorithm == "general_algorithm":
        filename_solution = "solution_part_b_n_" + str(n) + "_general.txt"
        filename_errors = "errors_n_" + str(n) + "_general.txt"
        filename_time = "timeused_" + str(n) + "_general.txt"
        system("rm" + " " + filename_solution + " " + filename_errors + " " + filename_time)
    if which_algorithm == "special_algorithm":
        filename_solution = "solution_part_b_n_" + str(n) + "_special.txt"
        filename_errors = "errors_n_" + str(n) + "_special.txt"
        filename_time = "timeused_" + str(n) + "_special.txt"
        system("rm" + " " + filename_solution + " " + filename_errors + " " + filename_time)
    if which_algorithm == "LU":
        filename_solution = "solution_part_b_n_" + str(n) + "_LU.txt"
        filename_errors = "errors_n_" + str(n) + "_LU.txt"
        filename_time = "timeused_" + str(n) + "_LU.txt"
        system("rm" + " " + filename_solution + " " + filename_errors + " " + filename_time)

print("Done")
