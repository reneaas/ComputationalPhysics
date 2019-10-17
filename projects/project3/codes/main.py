import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Library import StraightLine
from Library import PlottingTool
from statistics import mean
#plt.rc('text', usetex=True)

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

if compilation_instruction == "benchmark_laguerre":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")

    compile = str(input("Produce new data? Type yes or no: "))

    if compile == "yes":


        dimensions = [3,6]
        number_of_integration_points = [10*i for i in range(1,11)]

        #Runs the program for Gauss-Laguerre for both 3 and 6 dimensions, writing the results to file
        for i in dimensions:
            print("Running for %.f dimensions" % i)
            for n in number_of_integration_points:
                print("Executing for n = %.f" % n)
                outfilename = "dim_" + str(i) + "_" + str(n) + ".txt"
                integration_method = "2"
                os.system("./main.exe" + " " + outfilename + " " + integration_method + " " + str(i) + " " + str(n))



        integral_value = [[],[]]
        rel_error = [[],[]]
        timeused = [[],[]]

        #Reading the files generated above and read the data to lists.
        for i in range(len(dimensions)):
            for n in number_of_integration_points:
                filename = "dim_" + str(dimensions[i]) + "_" + str(n) + ".txt"
                with open(filename, "r") as infile:
                    line = infile.readline()
                    line = line.split()
                    integral_value[i].append(float(line[1]))
                    rel_error[i].append(float(line[2]))
                    timeused[i].append(float(line[3]))


                os.system("rm" + " " + filename)  #Delete .txt-files

        #Writing the combined data to main file for each dimension
        main_filename_3_dim = "benchmark_dim_3_.csv"
        main_filename_6_dim = "benchmark_dim_6_.csv"


        data_3_dim = {\
                "n" : number_of_integration_points,\
                "Integration value" : integral_value[0],\
                "Relative error" : rel_error[0], \
                "Time used" : timeused[0], \
                }

        data_6_dim = {\
                "n" : number_of_integration_points,\
                "Integration value" : integral_value[1],\
                "Relative error" : rel_error[1], \
                "Time used" : timeused[1], \
                }

        dataset3 = pd.DataFrame(data_3_dim)
        dataset6 = pd.DataFrame(data_6_dim)
        dataset3.to_csv(main_filename_3_dim, index = False)
        dataset6.to_csv(main_filename_6_dim, index = False)

        #Moves the files to appropriate file
        path = "results/laguerre";
        if not os.path.exists(path):
            os.makedirs(path)
        os.system("mv" + " " +  main_filename_3_dim + " " + main_filename_6_dim + " " + path)

    path = "results/laguerre/"

    #Reading the files generated above

    data_3_dim = pd.read_csv(path + "benchmark_dim_3_.csv", header = 0, names = ["n_3", "int_val_3", "rel_err_3", "time_3"])
    data_6_dim = pd.read_csv(path + "benchmark_dim_6_.csv", header = 0, names = ["n_6", "int_val_6", "rel_err_6", "time_6"])

    n = data_3_dim["n_3"]
    exact = 5*np.pi**2/(16*16)

    path = "results/laguerre/"

    data_3_dim = pd.read_csv(path + "benchmark_dim_3_.csv", header = 0, names = ["n_3", "int_val_3", "rel_err_3", "time_3"])
    data_6_dim = pd.read_csv(path + "benchmark_dim_6_.csv", header = 0, names = ["n_6", "int_val_6", "rel_err_6", "time_6"])

    n = data_3_dim["n_3"]
    exact = 5*np.pi**2/(16*16)

    #Comparing the results for 3 and 6 dimensions graphically
    figname1 = "integration_value.pdf"
    figname2 = "time_data.pdf"
    figname3 = "relative_error.pdf"


    plt.scatter(n, data_3_dim["int_val_3"], label="3 dimensions")
    plt.scatter(n, data_6_dim["int_val_6"], label="6 dimensions")
    plt.axhline(y = exact, ls = "--", label="Analytical value")
    plt.xlabel("n")
    plt.ylabel("Integration value")
    plt.legend()
    plt.savefig(figname1)
    plt.close()

    plt.plot(np.log10(n), np.log(data_3_dim["time_3"]), label="3 dimensions")
    plt.plot(np.log10(n), np.log(data_6_dim["time_6"]), label="6 dimensions")
    plt.xlabel("$log_{10}(n)$")
    plt.ylabel("$log_{10}(Time)$")
    plt.legend()
    plt.savefig(figname2)
    plt.close()


    plt.plot(n, (data_3_dim["rel_err_3"]), label="3 dimensions")
    plt.plot(n, (data_6_dim["rel_err_6"]), label="6 dimensions")
    plt.xlabel("n")
    plt.ylabel("$\epsilon_{relative}$")
    plt.legend()
    plt.savefig(figname3)
    plt.close()

    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " +  figname1 + " " + figname2 + " " + figname3 + " " + path)

if compilation_instruction == "compare_all":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    dimensions = 3
    number_of_integration_points = [5, 9, 15, 17, 21]

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


    number_of_monte_carlo_samples = [100, 1000, 10000, 100000, 1000000]
    #Next up is Brute force monte carlo, integration_method = "3".
    integration_method = "3"
    a = -2
    b = 2
    for n in number_of_monte_carlo_samples:
        print("Exectuting for Monte Carlo integration w/Brute force for n = " + str(n) + " samples")
        outfilename = "method_" + integration_method + "_" + str(n) + ".txt"
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b)
        os.system("./main.exe" + " " + arguments)


    integration_method = "4"
    max_radial_distance = 10
    for n in number_of_monte_carlo_samples:
        print("Exectuting for Monte Carlo integration w/importance sampling for n = " + str(n) + " samples")
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

    std_mean_MC_brute = []
    std_mean_MC_importance = []

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
    for n in number_of_monte_carlo_samples:
        infilename = "method_" + integration_method + "_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_MC_brute.append(float(numbers[0]))
            std_mean_MC_brute.append(float(numbers[1]))
            relative_error_MC_brute.append(float(numbers[2]))
            timeused_MC_brute.append(float(numbers[3]))
        os.system("rm" + " " + infilename)

    integration_method = "4"
    for n in number_of_monte_carlo_samples:
        infilename = "method_" + integration_method + "_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_MC_importance.append(float(numbers[0]))
            std_mean_MC_importance.append(float(numbers[1]))
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

    dataset_MC = {\
                    "N" : number_of_monte_carlo_samples,\
                    "$I_{bf}$" : integral_MC_brute,\
                    "$\sigma_{bf}$" : std_mean_MC_brute,\
                    "$I_{is}$" : integral_MC_importance,\
                    "$\sigma_{is}$" : std_mean_MC_importance\
                    }

    data_integrals = pd.DataFrame(dataset_integrals)
    data_relative_errors = pd.DataFrame(dataset_relative_errors)
    data_timeused = pd.DataFrame(dataset_timeused)
    data_MC = pd.DataFrame(dataset_MC)


    outfilename_integrals = "integrals_all_methods.txt"
    outfilename_relative_errors = "relative_error_all_methods.txt"
    outfilename_timeused = "timeused_all_methods.txt"
    outfilename_MC = "MC_comparison.txt"


    data_integrals.to_latex(outfilename_integrals, encoding='utf-8', escape = False, index = False)
    data_relative_errors.to_latex(outfilename_relative_errors, encoding='utf-8', escape = False, index = False)
    data_timeused.to_latex(outfilename_timeused, encoding='utf-8', escape = False, index = False)
    data_MC.to_latex(outfilename_MC, encoding='utf-8', escape = False, index = False)

    path = "results/benchmarks"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + outfilename_integrals + " " + outfilename_relative_errors + " " + outfilename_timeused + " " + outfilename_MC + " " + path)

if compilation_instruction == "ground_state":
    #Compiles the code without MPI.
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")

    number_of_monte_carlo_samples = [10**i for i in range(2,8)]
    integration_method = "4"
    Integrals = []
    relative_error = []



    for n in number_of_monte_carlo_samples:
        print("running code for n = " + str(n))
        outfilename = "monty_" + "n_" + str(n) + ".txt"
        os.system("./main.exe" + " " + outfilename + " " + integration_method + " " + str(n) )

    for n in number_of_monte_carlo_samples:
        infilename = "monty_" + "n_" + str(n) + ".txt"
        with open(infilename,"r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            Integrals.append(float(numbers[4]))
            relative_error.append(float(numbers[5]))
        os.system("rm" + " " + infilename)

    dataset = {"$N$" : number_of_monte_carlo_samples, "$\expval{H}$ [eV]" : Integrals, "$\epsilon$" : relative_error}

    dataset = pd.DataFrame(dataset)
    outfilename = "ground_state_energies.txt"
    path = "results/"
    dataset.to_latex(outfilename, encoding='utf-8', escape = False, index = False)

    if not os.path.exists(path):
        os.makedirs(path)

    os.system("mv" + " " + outfilename + " " + path)
