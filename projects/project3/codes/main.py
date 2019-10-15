import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Library import StraightLine
from Library import PlottingTool
plt.rc('text', usetex=True)

compilation_instruction = str(sys.argv[1])
if compilation_instruction == "no_mpi":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    print("executing")
    os.system("./main.exe")

if compilation_instruction == "0":
    print("executing")
    os.system("./main.exe")

if compilation_instruction == "mpi":
    print("compiling with mpi")
    os.system("mpicxx -O3 -c main_mpi.cpp")
    os.system("mpicxx -O3 -o main_mpi.exe main_mpi.o")
    print("executing with mpi")
    N = float(input("Specify number of monte carlo samples: "))
    os.system("mpirun -np 10 --oversubscribe ./main_mpi.exe" + " " + str(N))

if compilation_instruction == "mpi_timeit":
    print("Compiling main program WITH MPI...")
    os.system("mpicxx -O3 -c main_mpi.cpp")
    os.system("mpicxx -O3 -o main_mpi.exe main_mpi.o")
    #Number_of_monte_carlo_samples = [10**i for i in range(1,6)]
    Number_of_monte_carlo_samples = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]
    for N in Number_of_monte_carlo_samples:
        print("executing monte carlo integration WITH mpi for N = " + str(N) + " samples...")
        outfilename = "time_vs_n_montecarlo_mpi" + str(N) + ".txt"
        os.system("mpirun -np 10 --oversubscribe ./main_mpi.exe" + " " + str(N) + " " + "write_to_file" + " " + outfilename)

    #Lists to store the data.
    number_of_samples = []
    timeused_mpi = []
    computed_integrals_mpi = []
    timeused_no_mpi = []
    computed_integrals_no_mpi = []
    speedups = []
    relative_error_no_mpi = []
    relative_error_mpi = []


    for N in Number_of_monte_carlo_samples:
        infilename = "time_vs_n_montecarlo_mpi" + str(N) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            line = line.split()
            number_of_samples.append(float(line[0]))
            timeused_mpi.append(float(line[1]))
            computed_integrals_mpi.append(float(line[2]))
            relative_error_mpi.append(float(line[3]))


    #This part runs the same simulation but without mpi.
    print("Compiling main program WITHOUT MPI...")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    integration_method = "montecarlo_benchmarking"
    for N in Number_of_monte_carlo_samples:
        print("executing monte carlo integration WITHOUT mpi for N = " + str(N) + " samples...")
        outfilename = "time_vs_n_montecarlo_no_mpi_" + str(N) + ".txt"
        os.system("./main.exe" + " " + outfilename + " " + integration_method + " " + str(N))


    for N in Number_of_monte_carlo_samples:
        infilename = "time_vs_n_montecarlo_no_mpi_" + str(N) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            line = line.split()
            timeused_no_mpi.append(float(line[1]))
            computed_integrals_no_mpi.append(float(line[2]))
            relative_error_no_mpi.append(float(line[3]))

    #Computes speedup = T_1/T_p, where T_1 is time on a single processor, while T_p is the parallelized time.
    for T_1, T_p in zip(timeused_no_mpi, timeused_mpi):
        speedups.append(float(T_1/T_p))


    number_of_samples_log10 = [np.log10(N) for N in number_of_samples]            #Converts the number of samples N to log10(N) for readability.


    print("Writing the results to a file...")
    #Write all the results to congregate file.
    outfilename = "time_vs_n_montecarlo_mpi.txt"
    data = {\
            "$\log_{10}(N)$" : number_of_samples_log10,\
            "$\Delta t_\text{no MPI}$" : timeused_no_mpi,\
            "$\Delta t_\text{MPI}$" : timeused_mpi, \
            "$I_\text{no MPI}$" : computed_integrals_no_mpi, \
            "$I_\text{MPI}$" : computed_integrals_mpi, \
            "$T_1/T_p$" : speedups,\
            "Relative error, no mpi" : relative_error_no_mpi, \
            "Relative error, mpi" : relative_error_mpi \
            }
    dataset = pd.DataFrame(data)
    dataset.to_latex(outfilename, encoding='utf-8', escape = False, index = False)
    print("Results:")
    print("----------------------------------------------------------------------------------------------------------------------------------")
    print(dataset)
    print("----------------------------------------------------------------------------------------------------------------------------------")


    print("Moving the file to destination...")
    filename_table = "time_vs_n_montecarlo_mpi.txt"
    path = "results/monte_carlo/parallelization"                    #File destination
    if not os.path.exists(path):                                    #Creates directory if it doesn't exist.
        os.makedirs(path)
    os.system("mv" + " " + filename_table + " " + path)                #Moves the file to the correct directory.

    #Makes a plot of log10(N) vs speedup.

    figurename = "speedup_vs_n_montecarlo.pdf"
    xlabel = "$\log_{10}(N)$"
    ylabel = "$T_1/T_p$"
    labeltext = "Speedup = $T_1/T_p$"
    Line = StraightLine(x = number_of_samples_log10, y = speedups, number_of_datasets = 1)

    line = Line.straightline()
    Line.make_plot(labeltexts = labeltext, xlabel = xlabel, ylabel = ylabel, figurename = figurename)



    relative_error_mpi_log10 = [np.log10(i) for i in relative_error_mpi]
    relative_error_no_mpi_log10 = [np.log10(i) for i in relative_error_no_mpi]
    #Prepares the datasets to be plotted
    X_data = []
    Y_data = []
    X_data.append(number_of_samples_log10)
    X_data.append(number_of_samples_log10)
    Y_data.append(relative_error_no_mpi_log10)
    Y_data.append(relative_error_mpi_log10)

    #Creates an instance of PlottingTool and plots the datasets for relative error.
    Plotmaker = PlottingTool(x = X_data, y = Y_data, number_of_datasets = 2)
    figurename_relative_error = "relative_error_montecarlo.pdf"
    labeltexts = ["Without MPI", "With MPI"]
    xlabel = "$\log_{10} N$"
    ylabel = "$\log_{10} (\epsilon)$"
    type = "scatter"
    Plotmaker.plot(labeltexts = labeltexts, xlabel = xlabel, ylabel = ylabel, figurename = figurename_relative_error, type = type)

    figurename_relative_error_straightline = "relative_error_montecarlo_straightlines.pdf"
    Lines_relative_error = StraightLine(x = X_data, y = Y_data, number_of_datasets = 2)
    Lines_relative_error.straightline()
    Lines_relative_error.make_plot(labeltexts = labeltexts, xlabel = xlabel, ylabel = ylabel, figurename = figurename_relative_error_straightline)

    #Specifies filenames and labels for use with StraightLine.
    figurename_relative_error_with_mpi = "relative_error_mpi_montecarlo_straightline.pdf"
    labeltext_relative_error_mpi = "with MPI"
    figurename_relative_error_no_mpi = "relative_error_no_mpi_montecarlo_straightline.pdf"
    labeltext_relative_error_no_mpi = "without MPI"

    #Plots straight line and error bar for relative error with mpi
    Line_relative_error_mpi = StraightLine(x = number_of_samples_log10, y = relative_error_mpi_log10, number_of_datasets = 1)
    Line_relative_error_mpi.straightline()
    Line_relative_error_mpi.make_plot(labeltexts = labeltext_relative_error_mpi, xlabel = xlabel, ylabel = ylabel, figurename = figurename_relative_error_with_mpi)

    #Plots straigth line and error bar for relative error without mpi.
    Line_relative_error_no_mpi = StraightLine(x = number_of_samples_log10, y = relative_error_no_mpi_log10, number_of_datasets = 1)
    Line_relative_error_no_mpi.straightline()
    Line_relative_error_no_mpi.make_plot(labeltexts = labeltext_relative_error_no_mpi, xlabel = xlabel, ylabel = ylabel, figurename = figurename_relative_error_no_mpi)

    #Moves the file to destination
    path_figure = path
    os.system("mv" + " " + figurename + " " + figurename_relative_error + " " + figurename_relative_error_straightline + " " + figurename_relative_error_with_mpi + " " + figurename_relative_error_no_mpi + " " + path_figure)



    #Cleans up directory by removing excess files
    print("Cleaning up... Removing excess files...")
    for N in Number_of_monte_carlo_samples:
        filename_mpi = "time_vs_n_montecarlo_mpi" + str(N) + ".txt"
        filename_no_mpi = "time_vs_n_montecarlo_no_mpi_" + str(N) + ".txt"
        os.system("rm" + " " + filename_mpi + " " + filename_no_mpi)
    print("Done.")

if compilation_instruction == "benchmark_laguerre":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    dimensions = float(input("3 or 6 dimensions? "))
    number_of_integration_points = [5, 10, 15, 20, 25, 30, 35, 40]

    for n in number_of_integration_points:
        outfilename = "dim_" + str(dimensions) + "_" + str(n) + ".txt"
        integration_method = "2"
        os.system("./main.exe" + " " + outfilename + " " + integration_method + " " + str(dimensions) + " " + str(n))

    main_filename = "Benchmark_dim_" + str(dimensions) + "_.txt"

    integral_value = []
    rel_error = []
    timeused = []


    for n in number_of_integration_points:
        filename = "dim_" + str(dimensions) + "_" + str(n) + ".txt"
        with open(filename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            line = line.split()
            integral_value.append(float(line[0]))
            rel_error.append(float(line[2]))
            timeused.append(float(line[3]))


        os.system("rm" + " " + filename)  #Delete .txt-files

    data = {\
            "Integration value" : integral_value,\
            "Number of integration points" : number_of_integration_points,\
            "$\epsilon_\text{rel}$" : rel_error, \
            "Time used" : timeused, \
            }

    dataset = pd.DataFrame(data)
    dataset.to_latex(main_filename, encoding='utf-8', escape = False, index = False)

    print("Results:")
    print("----------------------------------------------------------------------------------------------------------------------------------")
    print(dataset)
    print("----------------------------------------------------------------------------------------------------------------------------------")

    path = "results/laguerre/benchmarks";
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " +  main_filename + " " + path)

if compilation_instruction == "compare_all":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    dimensions = 3
    number_of_integration_points = [5, 10, 15, 20, 25, 30]

    #First Gauss-Legendre; integration_method = "1".
    integration_method = "1"
    a = -3.0
    b = 3.0
    for n in number_of_integration_points:
        print("Exectuting for Gauss-Legendre method for n = " + str(n))
        outfilename = "method_" + integration_method + "_" + str(n) + ".txt"
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b)
        os.system("./main.exe" + " " + arguments)


    #Then Gauss-Laguerre, integration_method = "2"
    integration_method = "2"
    for n in number_of_integration_points:
        print("Exectuting for Gauss-Laguerre method for n = " + str(n))
        outfilename = "method_" + integration_method + "_" + str(n) + ".txt"
        os.system("./main.exe" + " " + outfilename + " " + integration_method + " " + str(dimensions) + " " + str(n))

    #Next up is Brute force monte carlo, integration_method = "3".
    integration_method = "3"
    for n in number_of_integration_points:
        print("Exectuting for Monte Carlo integration w/Brute force for n = " + str(n))
        outfilename = "method_" + integration_method + "_" + str(n) + ".txt"
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b)
        os.system("./main.exe" + " " + arguments)

    integration_method = "4"
    max_radial_distance = 10
    for n in number_of_integration_points:
        print("Exectuting for Monte Carlo integration w/importance sampling for n = " + str(n))
        outfilename = "method_" + integration_method + "_" + str(n) + ".txt"
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(max_radial_distance)
        os.system("./main.exe" + " " + arguments)



    timeused_gauleg = []
    timeused_gaulag = []
    timeused_MC_brute = []
    timeused_MC_importance = []

    relative_error_gauleg = []
    relative_error_gaulag = []
    relative_error_MC_brute = []
    relative_error_MC_importance = []

    integral_gauleg = []
    integral_gaulag = []
    integral_MC_brute = []
    integral_MC_importance = []

    integration_method = "1"
    for n in number_of_integration_points:
        infilename = "method_" + integration_method + "_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_gauleg.append(float(numbers[0]))
            relative_error_gauleg.append(float(numbers[2]))
            timeused_gauleg.append(float(numbers[3]))
        os.system("rm" + " " + infilename)

    integration_method = "2"
    for n in number_of_integration_points:
        infilename = "method_" + integration_method + "_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_gaulag.append(float(numbers[0]))
            relative_error_gaulag.append(float(numbers[2]))
            timeused_gaulag.append(float(numbers[3]))
        os.system("rm" + " " + infilename)

    integration_method = "3"
    for n in number_of_integration_points:
        infilename = "method_" + integration_method + "_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_MC_brute.append(float(numbers[0]))
            relative_error_MC_brute.append(float(numbers[2]))
            timeused_MC_brute.append(float(numbers[3]))
        os.system("rm" + " " + infilename)

    integration_method = "4"
    for n in number_of_integration_points:
        infilename = "method_" + integration_method + "_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_MC_importance.append(float(numbers[0]))
            relative_error_MC_importance.append(float(numbers[2]))
            timeused_MC_importance.append(float(numbers[3]))
        os.system("rm" + " " + infilename)


    #Create datasets.
    dataset_integrals = {\
                        "n" : number_of_integration_points, \
                        "Gauss-Legendre" : integral_gauleg, \
                        "Gauss-Laguerre" : integral_gaulag, \
                        "Monte Carlo Brute Force" : integral_MC_brute, \
                        "Monte carlo w/importance sampling" : integral_MC_importance \
                        }

    dataset_relative_errors = {\
                                "n" : number_of_integration_points, \
                                "Gauss-Legendre" : relative_error_gauleg, \
                                "Gauss-Laguerre" : relative_error_gaulag, \
                                "Monte Carlo w/Brute force" : relative_error_MC_brute, \
                                "Monte Carlo w/importance sampling" : relative_error_MC_importance\
                                }

    dataset_timeused = {\
                        "n" : number_of_integration_points, \
                        "Gauss-Legendre" : timeused_gauleg, \
                        "Gauss-Laguerre" : timeused_gaulag, \
                        "Monte Carlo w/Brute force" : timeused_MC_brute, \
                        "Monte Carlo w/importance sampling" : timeused_MC_importance\
                        }

    data_integrals = pd.DataFrame(dataset_integrals)
    data_relative_errors = pd.DataFrame(dataset_relative_errors)
    data_timeused = pd.DataFrame(dataset_timeused)


    outfilename_integrals = "integrals_all_methods.txt"
    outfilename_relative_errors = "relative_error_all_methods.txt"
    outfilename_timeused = "timeused_all_methods.txt"


    data_integrals.to_latex(outfilename_integrals, encoding='utf-8', escape = False, index = False)
    data_relative_errors.to_latex(outfilename_relative_errors, encoding='utf-8', escape = False, index = False)
    data_timeused.to_latex(outfilename_timeused, encoding='utf-8', escape = False, index = False)

    path = "results/benchmarks"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + outfilename_integrals + " " + outfilename_relative_errors + " " + outfilename_timeused + " " + path)

if compilation_instruction == "compare_montecarlo_methods":
    print("compiling code...")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    Number_of_monte_carlo_samples = [10, 100, 1000, 10000, 100000]

    #First run the code with brute force montecarlo, integration method "3":
    integration_method = "3"
    n = 10
    a = -10
    b = 10
    for N in Number_of_monte_carlo_samples:
        outfilename = "bruteforce_N_" + str(N) + ".txt"
        os.system("./main.exe" + " " + outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b) + " " + str(N))

    #Then we run the code for monte carlo with importance sampling, integration_method = "4"
    integration_method = "4"
    n = 10
    max_radial_distance = 10;
    for N in Number_of_monte_carlo_samples:
        outfilename = "importancesampling_N_" + str(N) + ".txt"
        os.system("./main.exe" + " " + outfilename + " " + integration_method + " " + str(n) + " " + str(max_radial_distance) + " " + str(N))

    #Defines empty lists to store data
    integrals_bruteforce = []
    timeused_bruteforce = []
    relative_error_bruteforce = []
    standard_deviation_bruteforce = []

    integrals_ImportanceSampling = []
    timeused_ImportanceSampling = []
    relative_error_ImportanceSampling = []
    standard_deviation_ImportanceSampling = []

    #Read the files with the results for brute force monte carlo integration
    for N in Number_of_monte_carlo_samples:
        infilename = "bruteforce_N_" + str(N) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            values = line.split()
            integrals_bruteforce.append(float(values[0]))
            relative_error_bruteforce.append(float(values[2]))
            timeused_bruteforce.append(float(values[3]))
