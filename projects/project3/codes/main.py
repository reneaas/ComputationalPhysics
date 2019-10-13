import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
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
    Number_of_monte_carlo_samples = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]
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
    plt.scatter(number_of_samples_log10, speedups, label = r"Speedup = $T_1/T_p$")
    plt.xlabel(r"$\log_{10}(N)$")
    plt.ylabel(r"$T_1/T_p$")
    plt.legend()
    plt.savefig(figurename)
    plt.close()

    #Moves the file to destination
    path_figure = path
    os.system("mv" + " " + figurename + " " + path_figure)



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
