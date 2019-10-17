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

    number_of_monte_carlo_samples = [10**i for i in range(2,8)]
    #number_of_monte_carlo_samples = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 100000000]
    max_radial_distance = 10
    a = -3
    b = 3


    #Runs the code with MPI in Cartesian coordinates.
    integration_method = "1"
    for n in number_of_monte_carlo_samples:
        print("Cartesian coordinates WITH MPI for n = " + str(n) + " samples...")
        outfilename = "MPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b)
        os.system("mpirun -np 2 --oversubscribe ./main_mpi.exe" + " " + arguments)


    #Runs the code with MPI in spherical coordinates.
    integration_method = "2"
    for n in number_of_monte_carlo_samples:
        print("Spherical coordinates WITH mpi for n = " + str(n) + " samples...")
        outfilename = "MPI_integrationmethod_" + integration_method + "_spherical_n_" + str(n) + ".txt"
        arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(max_radial_distance)
        os.system("mpirun -np 2 --oversubscribe ./main_mpi.exe" + " " + arguments)

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
                    "MPI (BF)" : integral_MPI_cartesian,\
                    "MPI (IS)" : integral_MPI_spherical,\
                    "No MPI (BF)" : integral_cartesian,\
                    "No MPI (IS)" : integral_spherical\
                }

    STDmean = {\
                    "n" : number_of_monte_carlo_samples,\
                    "$\sigma_\text{BF}$ (MPI)" : STDmean_MPI_cartesian,\
                    "$\sigma_\text{IS}$ (MPI)" : STDmean_MPI_spherical,\
                    "$\sigma_\text{BF}$" : STDmean_cartesian,\
                    "$\sigma_\text{IS}$" : STDmean_spherical\
                }

    RelativeError = {\
                    "n" : number_of_monte_carlo_samples,\
                    "$\epsilon_\text{BF}$ (MPI)" : RelativeError_MPI_cartesian,\
                    "$\epsilon_\text{IS}$ (MPI)" : RelativeError_MPI_spherical,\
                    "$\epsilon_\text{BF}$" : RelativeError_cartesian,\
                    "$\epsilon_\text{IS}$ (IS)" : RelativeError_spherical\
                }

    timeused = {\
                    "n" : number_of_monte_carlo_samples,\
                    "$\Delta t_\text{BF}$ (MPI)" : timeused_MPI_cartesian,\
                    "$\Delta t_\text{IS}$ (MPI)" : timeused_MPI_spherical,\
                    "$\Delta t_\text{BF}$" : timeused_cartesian,\
                    "$\Delta t_\text{IS}$ (IS)" : timeused_spherical,\
                }

    speedup = {\
                    "n" : number_of_monte_carlo_samples,\
                    "speedup (BF)" : speedup_cartesian,\
                    "speedup (IS)" : speedup_spherical,\
                }

    integrals = pd.DataFrame(integrals)
    STDmean = pd.DataFrame(STDmean)
    RelativeError = pd.DataFrame(RelativeError)
    timeused = pd.DataFrame(timeused)
    speedup = pd.DataFrame(speedup)

    outfilename_integrals = "integrals_MCmethods.txt"
    outfilename_STDmean = "STDmean_MCmethods.txt"
    outfilename_RelativeError = "RelativeError_MCmethods.txt"
    outfilename_timeused = "timeused_MCmethods.txt"
    outfilename_speedup = "speedup_MCmethods.txt"

    integrals.to_latex(outfilename_integrals, encoding='utf-8', escape = False, index = False)
    STDmean.to_latex(outfilename_STDmean, encoding='utf-8', escape = False, index = False)
    RelativeError.to_latex(outfilename_RelativeError, encoding='utf-8', escape = False, index = False)
    timeused.to_latex(outfilename_timeused, encoding='utf-8', escape = False, index = False)
    speedup.to_latex(outfilename_speedup, encoding='utf-8', escape = False, index = False)

    path = "results/monte_carlo/benchmarks"
    if not os.path.exists(path):
        os.makedirs(path)


    os.system("mv" + " " + outfilename_integrals + " " + path)
    os.system("mv" + " " + outfilename_STDmean + " " + path)
    os.system("mv" + " " + outfilename_RelativeError + " " + path)
    os.system("mv" + " " + outfilename_timeused + " " + path)
    os.system("mv" + " " + outfilename_speedup + " " + path)


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


    X_data = []
    for i in range(4):
        X_data.append(np.log10(number_of_monte_carlo_samples))

    timeused_mat = []
    timeused_mat.append(np.log10(timeused_MPI_cartesian))
    timeused_mat.append(np.log10(timeused_MPI_spherical))
    timeused_mat.append(np.log10(timeused_cartesian))
    timeused_mat.append(np.log10(timeused_spherical))

    figurename_lines_timeused = "Time_used_for_all_Monte_Carlo_methods.pdf"
    labeltexts_lines = ["Brute force w/MPI", "Importance sampling w/MPI", "Brute force", "Importance sampling"]
    xlabel_lines = "$\log_{10}(N)$"
    ylabel_lines_timeused = "$\log_{10}(t)$"
    lines_timeused = StraightLine(X_data, timeused_mat, 4)
    lines_timeused.straightline()
    lines_timeused.make_plot(labeltexts_lines, xlabel_lines, ylabel_lines_timeused, figurename_lines_timeused)



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

if compilation_instruction == "multiple_MC":



    print("Compiling main program WITH MPI...")
    os.system("mpicxx -O3 -c main_mpi.cpp")
    os.system("mpicxx -O3 -o main_mpi.exe main_mpi.o")

    #Compiles the code without MPI.
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")

    number_of_monte_carlo_samples = [10**i for i in range(2,8)]
    #number_of_monte_carlo_samples = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 100000000]
    max_radial_distance = 10
    a = -2
    b = 2


    #How many times we want to run the code for the same value of monte carlo sample
    m_runs = 100

    #Dictionaries to hold calculated values
    dict_Integral_BF = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Integral_BF_MPI = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Integral_IS = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Integral_IS_MPI = {str(i):[] for i in number_of_monte_carlo_samples}

    dict_Error_BF = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Error_BF_MPI = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Error_IS = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Error_IS_MPI = {str(i):[] for i in number_of_monte_carlo_samples}

    dict_Std_BF = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Std_BF_MPI = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Std_IS = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Std_IS_MPI = {str(i):[] for i in number_of_monte_carlo_samples}

    dict_Time_BF = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Time_BF_MPI = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Time_IS = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Time_IS_MPI = {str(i):[] for i in number_of_monte_carlo_samples}

    dict_Speedup_BF = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Speedup_IS = {str(i):[] for i in number_of_monte_carlo_samples}

    dict_Ground_state_IS = {str(i):[] for i in number_of_monte_carlo_samples}
    dict_Error_Ground_state_IS = {str(i):[] for i in number_of_monte_carlo_samples}

    #Multiple runs of same montecarlo value
    for m in range(0,m_runs):

        #Runs the code with MPI in Cartesian coordinates.
        integration_method = "1"
        for n in number_of_monte_carlo_samples:
            outfilename = "MPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
            arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b)
            os.system("mpirun -np 2 --oversubscribe ./main_mpi.exe" + " " + arguments)

        #Adding values to dictionary
        for n in number_of_monte_carlo_samples:
            infilename = "MPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
            with open(infilename,"r") as infile:
                lines = infile.readlines()
                line = lines[0]
                numbers = line.split()
                dict_Integral_BF_MPI[str(n)].append(float(numbers[0]))
                dict_Std_BF_MPI[str(n)].append(float(numbers[1]))
                dict_Error_BF_MPI[str(n)].append(float(numbers[2]))
                dict_Time_BF_MPI[str(n)].append(float(numbers[3]))
            os.system("rm" + " " + infilename)

        #Runs the code with MPI in spherical coordinates.
        integration_method = "2"
        for n in number_of_monte_carlo_samples:
            outfilename = "MPI_integrationmethod_" + integration_method + "_spherical_n_" + str(n) + ".txt"
            arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(max_radial_distance)
            os.system("mpirun -np 2 --oversubscribe ./main_mpi.exe" + " " + arguments)

        #Adding values to dictionary
        for n in number_of_monte_carlo_samples:
            infilename = "MPI_integrationmethod_" + integration_method + "_spherical_n_" + str(n) + ".txt"
            with open(infilename,"r") as infile:
                lines = infile.readlines()
                line = lines[0]
                numbers = line.split()
                dict_Integral_IS_MPI[str(n)].append(float(numbers[0]))
                dict_Std_IS_MPI[str(n)].append(float(numbers[1]))
                dict_Error_IS_MPI[str(n)].append(float(numbers[2]))
                dict_Time_IS_MPI[str(n)].append(float(numbers[3]))
            os.system("rm" + " " + infilename)

        #Runs the code without MPI in Cartesian coordinates
        integration_method = "3"
        for n in number_of_monte_carlo_samples:
            outfilename = "NoMPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
            arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b)
            os.system("./main.exe" + " " + arguments)

        #Adding values to dictionary
        for n in number_of_monte_carlo_samples:
            infilename = outfilename = "NoMPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
            with open(infilename,"r") as infile:
                lines = infile.readlines()
                line = lines[0]
                numbers = line.split()
                dict_Integral_BF[str(n)].append(float(numbers[0]))
                dict_Std_BF[str(n)].append(float(numbers[1]))
                dict_Error_BF[str(n)].append(float(numbers[2]))
                dict_Time_BF[str(n)].append(float(numbers[3]))
            os.system("rm" + " " + infilename)

        #Runs the code without MPI in spherical coordinates
        integration_method = "4"
        for n in number_of_monte_carlo_samples:
            outfilename = "NoMPI_integrationmethod_" + integration_method + "_spherical_n_" + str(n) + ".txt"
            arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(max_radial_distance)
            os.system("./main.exe" + " " + arguments)

        #Adding values to dictionary
        for n in number_of_monte_carlo_samples:
            infilename = "NoMPI_integrationmethod_" + integration_method + "_spherical_n_" + str(n) + ".txt"
            with open(infilename,"r") as infile:
                lines = infile.readlines()
                line = lines[0]
                numbers = line.split()
                dict_Integral_IS[str(n)].append(float(numbers[0]))
                dict_Std_IS[str(n)].append(float(numbers[1]))
                dict_Error_IS[str(n)].append(float(numbers[2]))
                dict_Time_IS[str(n)].append(float(numbers[3]))
                dict_Ground_state_IS[str(n)].append(float(numbers[4]))
                dict_Error_Ground_state_IS[str(n)].append(float(numbers[5]))
            os.system("rm" + " " + infilename)
    print("--------------------------------------------------------------------------------------------------------")
    print(dict_Integral_BF)
    print("--------------------------------------------------------------------------------------------------------")


    for n in number_of_monte_carlo_samples:
        key = str(n)
        dict_Integral_BF[key] = np.mean(dict_Integral_BF[key])
        dict_Integral_BF_MPI[key] = np.mean(dict_Integral_BF_MPI[key])
        dict_Integral_IS[key] = np.mean(dict_Integral_IS[key])
        dict_Integral_IS_MPI[key] = np.mean(dict_Integral_IS_MPI[key])

        dict_Error_BF[key] = np.mean(dict_Error_BF[key])
        dict_Error_BF_MPI[key] = np.mean(dict_Error_BF_MPI[key])
        dict_Error_IS[key] = np.mean(dict_Error_IS[key])
        dict_Error_IS_MPI[key] = np.mean(dict_Error_IS_MPI[key])

        dict_Std_BF[key] = np.mean(dict_Std_BF[key])
        dict_Std_BF_MPI[key] = np.mean(dict_Std_BF_MPI[key])
        dict_Std_IS[key] = np.mean(dict_Std_IS[key])
        dict_Std_IS_MPI[key] = np.mean(dict_Std_IS_MPI[key])

        dict_Time_BF[key] = np.mean(dict_Time_BF[key])
        dict_Time_BF_MPI[key] = np.mean(dict_Time_BF_MPI[key])
        dict_Time_IS[key] = np.mean(dict_Time_IS[key])
        dict_Time_IS_MPI[key] = np.mean(dict_Time_IS_MPI[key])


        dict_Ground_state_IS[key] = np.mean(dict_Ground_state_IS[key])
        dict_Error_Ground_state_IS[key] = np.mean(dict_Error_Ground_state_IS[key])

        dict_Speedup_BF[key] = float(dict_Time_BF[key])/float(dict_Time_BF_MPI[key])
        dict_Speedup_IS[key] = float(dict_Time_IS[key])/float(dict_Time_IS_MPI[key])


    I_BF = [dict_Integral_BF[str(i)] for i in number_of_monte_carlo_samples]
    I_BF_MPI = [dict_Integral_BF_MPI[str(i)] for i in number_of_monte_carlo_samples]
    I_IS = [dict_Integral_IS[str(i)] for i in number_of_monte_carlo_samples]
    I_IS_MPI = [dict_Integral_IS_MPI[str(i)] for i in number_of_monte_carlo_samples]

    E_BF = [dict_Error_BF[str(i)] for i in number_of_monte_carlo_samples]
    E_BF_MPI = [dict_Error_BF_MPI[str(i)] for i in number_of_monte_carlo_samples]
    E_IS = [dict_Error_IS[str(i)] for i in number_of_monte_carlo_samples]
    E_IS_MPI = [dict_Error_BF_MPI[str(i)] for i in number_of_monte_carlo_samples]

    Std_BF = [dict_Std_BF[str(i)] for i in number_of_monte_carlo_samples]
    Std_BF_MPI = [dict_Std_BF_MPI[str(i)] for i in number_of_monte_carlo_samples]
    Std_IS = [dict_Std_IS[str(i)] for i in number_of_monte_carlo_samples]
    Std_IS_MPI = [dict_Std_IS_MPI[str(i)] for i in number_of_monte_carlo_samples]

    Time_BF = [dict_Time_BF[str(i)] for i in number_of_monte_carlo_samples]
    Time_BF_MPI = [dict_Time_BF_MPI[str(i)] for i in number_of_monte_carlo_samples]
    Time_IS = [dict_Time_IS[str(i)] for i in number_of_monte_carlo_samples]
    Time_IS_MPI = [dict_Time_IS_MPI[str(i)] for i in number_of_monte_carlo_samples]

    Ground_state_IS = [dict_Ground_state_IS[str(i)] for i in number_of_monte_carlo_samples]
    E_Ground_state_IS = [dict_Error_Ground_state_IS[str(i)] for i in number_of_monte_carlo_samples]

    Speedup_BF = [dict_Speedup_BF[str(i)] for i in number_of_monte_carlo_samples]
    Speedup_IS = [dict_Speedup_IS[str(i)] for i in number_of_monte_carlo_samples]



    Integrals = {\
                    "$N$" : number_of_monte_carlo_samples,\
                    "$I_\text{BF}$" : I_BF,\
                    "$I_\text{BF}$ (MPI)" : I_BF_MPI,\
                    "$I_\text{IS}$" : I_IS,\
                    "$I_\text{IS}$ (MPI)$" : I_IS_MPI \
                }

    Errors = {\
                    "$N$" : number_of_monte_carlo_samples,\
                    "$\epsilon_\text{BF}$" : E_BF,\
                    "$\epsilon_\text{BF}$ (MPI)" : E_BF_MPI,\
                    "$\epsilon_\text{IS}$" : E_IS,\
                    "$\epsilon_\text{IS}$ (MPI)$" : E_IS_MPI\
                }

    Standard_deviations = {\
                    "$N$" : number_of_monte_carlo_samples,\
                    "$\sigma_\text{BF}$" : Std_BF,\
                    "$\sigma_\text{BF}$ (MPI)" : Std_BF_MPI,\
                    "$\sigma_\text{IS}$" : Std_IS,\
                    "$\sigma_\text{IS}$ (MPI)$" : Std_IS_MPI\
                }

    Times = {\
                    "$N$" : number_of_monte_carlo_samples,\
                    "$t_\text{BF}$" : Time_BF,\
                    "$t_\text{BF}$ (MPI)" : Time_BF_MPI,\
                    "$t_\text{IS}$" : Time_IS,\
                    "$t_\text{IS}$ (MPI)$" : Time_IS_MPI\
                }

    Ground_state = {\
                    "$N$" : number_of_monte_carlo_samples,\
                    "$\expval{H}$" : Ground_state_IS,\
                    "$\epsilon$" : E_Ground_state_IS\
                }

    Speedup = {\
               "$N$" : number_of_monte_carlo_samples,\
               "$T_1/T_p$ (BF)" : Speedup_BF,\
               "$T_1/T_p$ (IS)" : Speedup_IS\

             }

    Integrals = pd.DataFrame(Integrals)
    Errors = pd.DataFrame(Errors)
    Standard_deviations = pd.DataFrame(Standard_deviations)
    Times = pd.DataFrame(Times)
    Ground_state = pd.DataFrame(Ground_state)
    Speedup = pd.DataFrame(Speedup)


    outfilename_Integrals = "MC_Integrals.txt"
    outfilename_Errors = "MC_Errors.txt"
    outfilename_Standard_deviations = "MC_Standard_deviations.txt"
    outfilename_Times = "MC_Times.txt"
    outfilename_Ground_state = "MC_Ground_state.txt"
    outfilename_Speedup = "MC_Speedup.txt"

    path = "results/monte_carlo/benchmarks"
    if not os.path.exists(path):
        os.makedirs(path)


    Integrals.to_latex(outfilename_Integrals, encoding='utf-8', escape = False, index = False)
    Errors.to_latex(outfilename_Errors, encoding='utf-8', escape = False, index = False)
    Standard_deviations.to_latex(outfilename_Standard_deviations, encoding='utf-8', escape = False, index = False)
    Times.to_latex(outfilename_Times, encoding='utf-8', escape = False, index = False)
    Ground_state.to_latex(outfilename_Ground_state, encoding='utf-8', escape = False, index = False)
    Speedup.to_latex(outfilename_Speedup, encoding='utf-8', escape = False, index = False)

    print("Integrals:")
    print("------------------------------------------------------------------------------------")
    print(Integrals)
    print("------------------------------------------------------------------------------------")
    print("Ground state")
    print("------------------------------------------------------------------------------------")
    print(Ground_state)
    print("------------------------------------------------------------------------------------")
