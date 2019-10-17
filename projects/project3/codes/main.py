import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Library import StraightLine
from Library import PlottingTool
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
    dimensions = [3,6]
    number_of_integration_points = [5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50]


    for i in dimensions:
        for n in number_of_integration_points:
            print("Executing for n = %.f" % n)
            outfilename = "dim_" + str(i) + "_" + str(n) + ".txt"
            integration_method = "2"
            os.system("./main.exe" + " " + outfilename + " " + integration_method + " " + str(i) + " " + str(n))



    integral_value = [[],[]]
    rel_error = [[],[]]
    timeused = [[],[]]


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


    path = "results/laguerre";
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " +  main_filename_3_dim + " " + main_filename_6_dim + " " + path)

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

if compilation_instruction == "compare_MC":
    #Compiles the code for with MPI
    print("Compiling main program WITH MPI...")
    os.system("mpicxx -O3 -c main_mpi.cpp")
    os.system("mpicxx -O3 -o main_mpi.exe main_mpi.o")

    number_of_monte_carlo_samples = [10**i for i in range(1,9)]
    max_radial_distance = 10
    a = -3
    b = 3


    #Runs the code with MPI in Cartesian coordinates.
    integration_method = "1"
    for n in number_of_monte_carlo_samples:
        print("Cartesian coordinates WITH MPI for n = " + str(n) + " samples...")
        outfilename = "MPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b)
        os.system("mpirun -np 10 --oversubscribe ./main_mpi.exe" + " " + arguments)


    #Runs the code with MPI in spherical coordinates.
    integration_method = "2"
    for n in number_of_monte_carlo_samples:
        print("Spherical coordinates WITH mpi for n = " + str(n) + " samples...")
        outfilename = "MPI_integrationmethod_" + integration_method + "_spherical_n_" + str(n) + ".txt"
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(max_radial_distance)
        os.system("mpirun -np 10 --oversubscribe ./main_mpi.exe" + " " + arguments)

    #Compiles the code without MPI.
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")

    #Runs the code without MPI in Cartesian coordinates
    integration_method = "3"
    for n in number_of_monte_carlo_samples:
        outfilename = "NoMPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
        print("running code without MPI for n = " + str(n) + " samples ...")
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b)
        os.system("./main.exe" + " " + arguments)


    #Runs the code without MPI in spherical coordinates
    integration_method = "4"
    for n in number_of_monte_carlo_samples:
        outfilename = "NoMPI_integrationmethod_" + integration_method + "_spherical_n_" + str(n) + ".txt"
        print("running code without MPI for n = " + str(n) + " samples ...")
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(max_radial_distance)
        os.system("./main.exe" + " " + arguments)


    integral_MPI_spherical = []
    STDmean_MPI_spherical = []
    RelativeError_MPI_spherical = []
    timeused_MPI_spherical = []

    integral_MPI_cartesian = []
    STDmean_MPI_cartesian = []
    RelativeError_MPI_cartesian = []
    timeused_MPI_cartesian = []

    integral_spherical = []
    STDmean_spherical = []
    RelativeError_spherical = []
    timeused_spherical = []

    integral_cartesian = []
    STDmean_cartesian = []
    RelativeError_cartesian = []
    timeused_cartesian = []

    speedup_spherical = []
    speedup_cartesian = []


    #Read the computed data for MPI in cartesian coordinates
    integration_method = "1"
    for n in number_of_monte_carlo_samples:
        infilename = "MPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_MPI_cartesian.append(float(numbers[0]))
            STDmean_MPI_cartesian.append(float(numbers[1]))
            RelativeError_MPI_cartesian.append(float(numbers[2]))
            timeused_MPI_cartesian.append(float(numbers[3]))
        os.system("rm" + " " + infilename)

    #Read the computed data for MPI in spherical coordinates
    integration_method = "2"
    for n in number_of_monte_carlo_samples:
        infilename = "MPI_integrationmethod_" + integration_method + "_spherical_n_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_MPI_spherical.append(float(numbers[0]))
            STDmean_MPI_spherical.append(float(numbers[1]))
            RelativeError_MPI_spherical.append(float(numbers[2]))
            timeused_MPI_spherical.append(float(numbers[3]))
        os.system("rm" + " " + infilename)


    #Read the computed data without MPI in Cartesian coordinates
    integration_method = "3"
    for n in number_of_monte_carlo_samples:
        infilename = "NoMPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_cartesian.append(float(numbers[0]))
            STDmean_cartesian.append(float(numbers[1]))
            RelativeError_cartesian.append(float(numbers[2]))
            timeused_cartesian.append(float(numbers[3]))
        os.system("rm" + " " + infilename)


    #Read the computed data without MPI in spherical coordinates
    integration_method = "4"
    for n in number_of_monte_carlo_samples:
        infilename = "NoMPI_integrationmethod_" + integration_method + "_spherical_n_" + str(n) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            line = lines[0]
            numbers = line.split()
            integral_spherical.append(float(numbers[0]))
            STDmean_spherical.append(float(numbers[1]))
            RelativeError_spherical.append(float(numbers[2]))
            timeused_spherical.append(float(numbers[3]))
        os.system("rm" + " " + infilename)


    for time_mpi_spherical, time_nompi_spherical in zip(timeused_MPI_spherical, timeused_spherical):
        speedup_spherical.append(time_nompi_spherical/time_mpi_spherical)

    for time_mpi_cartesian, time_nompi_cartesian in zip(timeused_MPI_cartesian, timeused_cartesian):
        speedup_cartesian.append(time_nompi_cartesian/time_mpi_cartesian)

    integrals = {\
                    "n" : number_of_monte_carlo_samples,\
                    "MPI (cartesian)" : integral_MPI_cartesian,\
                    "MPI (spherical)" : integral_MPI_spherical,\
                    "No MPI (cartesian)" : integral_cartesian,\
                    "No MPI (spherical)" : integral_spherical\
                }

    STDmean = {\
                    "n" : number_of_monte_carlo_samples,\
                    "MPI (cartesian)" : STDmean_MPI_cartesian,\
                    "MPI (spherical)" : STDmean_MPI_spherical,\
                    "No MPI (cartesian)" : STDmean_cartesian,\
                    "No MPI (spherical)" : STDmean_spherical\
                }

    RelativeError = {\
                    "n" : number_of_monte_carlo_samples,\
                    "MPI (cartesian)" : RelativeError_MPI_cartesian,\
                    "MPI (spherical)" : RelativeError_MPI_spherical,\
                    "No MPI (cartesian)" : RelativeError_cartesian,\
                    "No MPI (spherical)" : RelativeError_spherical\
                }

    timeused = {\
                    "n" : number_of_monte_carlo_samples,\
                    "MPI (cartesian)" : timeused_MPI_cartesian,\
                    "MPI (spherical)" : timeused_MPI_spherical,\
                    "No MPI (cartesian)" : timeused_cartesian,\
                    "No MPI (spherical)" : timeused_spherical,\
                }

    speedup = {\
                    "n" : number_of_monte_carlo_samples,\
                    "speedup (cartesian)" : speedup_cartesian,\
                    "speedup (spherical)" : speedup_spherical,\
                }

    integrals = pd.DataFrame(integrals)
    STDmean = pd.DataFrame(STDmean)
    RelativeError = pd.DataFrame(RelativeError)
    timeused = pd.DataFrame(timeused)
    speedup = pd.DataFrame(speedup)

    print("----------------------------------------------------------------------------------")
    print("integrals:")
    print("----------------------------------------------------------------------------------")
    print(integrals)
    print("----------------------------------------------------------------------------------")
    print("STDmean:")
    print("----------------------------------------------------------------------------------")
    print(STDmean)
    print("----------------------------------------------------------------------------------")
    print("Relative error:")
    print("----------------------------------------------------------------------------------")
    print(RelativeError)
    print("----------------------------------------------------------------------------------")
    print("Time used:")
    print("----------------------------------------------------------------------------------")
    print(timeused)
    print("----------------------------------------------------------------------------------")
    print("Speedup:")
    print("----------------------------------------------------------------------------------")
    print(speedup)


if compilation_instruction == "compare_laguerre":

    path = "results/laguerre/"

    data_3_dim = pd.read_csv(path + "benchmark_dim_3_.csv", header = 0, names = ["n_3", "int_val_3", "rel_err_3", "time_3"])
    data_6_dim = pd.read_csv(path + "benchmark_dim_6_.csv", header = 0, names = ["n_6", "int_val_6", "rel_err_6", "time_6"])

    n = data_3_dim["n_3"]
    exact = 5*np.pi**2/(16*16)



    plt.scatter(n, data_3_dim["int_val_3"], label="3 dimensions")
    plt.scatter(n, data_6_dim["int_val_6"], label="6 dimensions")
    plt.axhline(y = exact, ls = "--", label="Analytical value")
    plt.legend()
    plt.show()


    plt.plot(n, np.log(data_3_dim["time_3"]), label="3 dimensions")
    plt.plot(n, np.log(data_6_dim["time_6"]), label="6 dimensions")
    plt.legend()
    plt.show()


    plt.plot(n, (data_3_dim["rel_err_3"]), label="3 dimensions")
    plt.plot(n, (data_6_dim["rel_err_6"]), label="6 dimensions")
    plt.legend()
    plt.show()
