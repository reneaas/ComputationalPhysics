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

if compilation_instruction == "multiple_MC":

    compile = str(input("Produce new data? Type yes or no: "))

    if compile == "yes":

        m_runs = int(input("How many sets of data do you want?: "))

        print("Compiling main program WITH MPI...")
        os.system("mpicxx -O3 -c main_mpi.cpp")
        os.system("mpicxx -O3 -o main_mpi.exe main_mpi.o lib.o")

        #Compiles the code without MPI.
        print("compiling main progam without MPI....")
        os.system("c++ -O3 -c main.cpp lib.cpp")
        os.system("c++ -O3 -o main.exe main.cpp lib.o")

        number_of_monte_carlo_samples = [100*i for i in range(1,101)]
        #number_of_monte_carlo_samples = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 100000000]
        max_radial_distance = 10
        a = -2.7
        b = 2.7


        #How many times we want to run the code for the same value of monte carlo samples
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
        print("executing...")
        for m in range(0,m_runs):
            print("Iteration number =", m+1, " of", m_runs)

            #Runs the code with MPI in Cartesian coordinates.
            integration_method = "1"
            for n in number_of_monte_carlo_samples:
                outfilename = "MPI_integrationmethod_" + integration_method + "_cartesian_n_" + str(n) + ".txt"
                arguments = outfilename + " " + integration_method + " " + str(n) + " " + str(a) + " " + str(b)
                os.system("mpirun -np 4 --oversubscribe ./main_mpi.exe" + " " + arguments)

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
                os.system("mpirun -np 4 --oversubscribe ./main_mpi.exe" + " " + arguments)

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
        E_IS_MPI = [dict_Error_IS_MPI[str(i)] for i in number_of_monte_carlo_samples]

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
                        "N" : number_of_monte_carlo_samples,\
                        "$BF" : I_BF,\
                        "BF (MPI)" : I_BF_MPI,\
                        "IS" : I_IS,\
                        "IS (MPI)$" : I_IS_MPI \
                    }

        Errors = {\
                        "N" : number_of_monte_carlo_samples,\
                        "BF" : E_BF,\
                        "BF (MPI)" : E_BF_MPI,\
                        "IS" : E_IS,\
                        "IS (MPI)" : E_IS_MPI\
                    }

        Standard_deviations = {\
                        "N" : number_of_monte_carlo_samples,\
                        "BF" : Std_BF,\
                        "BF (MPI)" : Std_BF_MPI,\
                        "IS" : Std_IS,\
                        "IS (MPI)" : Std_IS_MPI\
                    }

        Times = {\
                        "N" : number_of_monte_carlo_samples,\
                        "BF" : Time_BF,\
                        "BF (MPI)" : Time_BF_MPI,\
                        "IS" : Time_IS,\
                        "IS (MPI)" : Time_IS_MPI\
                    }

        Ground_state = {\
                        "N" : number_of_monte_carlo_samples,\
                        "H" : Ground_state_IS,\
                        "Error" : E_Ground_state_IS\
                    }

        Speedup = {\
                   "N" : number_of_monte_carlo_samples,\
                   "Speedup (BF)" : Speedup_BF,\
                   "Speedup (IS)" : Speedup_IS\

                 }

        Integrals = pd.DataFrame(Integrals)
        Errors = pd.DataFrame(Errors)
        Standard_deviations = pd.DataFrame(Standard_deviations)
        Times = pd.DataFrame(Times)
        Ground_state = pd.DataFrame(Ground_state)
        Speedup = pd.DataFrame(Speedup)


        outfilename_Integrals = "MC_Integrals.csv"
        outfilename_Errors = "MC_Errors.csv"
        outfilename_Standard_deviations = "MC_Standard_deviations.csv"
        outfilename_Times = "MC_Times.csv"
        outfilename_Ground_state = "MC_Ground_state.csv"
        outfilename_Speedup = "MC_Speedup.csv"

        """
        Integrals.to_latex(outfilename_Integrals, encoding='utf-8', escape = False, index = False)
        Errors.to_latex(outfilename_Errors, encoding='utf-8', escape = False, index = False)
        Standard_deviations.to_latex(outfilename_Standard_deviations, encoding='utf-8', escape = False, index = False)
        Times.to_latex(outfilename_Times, encoding='utf-8', escape = False, index = False)
        Ground_state.to_latex(outfilename_Ground_state, encoding='utf-8', escape = False, index = False)
        Speedup.to_latex(outfilename_Speedup, encoding='utf-8', escape = False, index = False)
        """

        Integrals.to_csv(outfilename_Integrals, index = False)
        Errors.to_csv(outfilename_Errors, index = False)
        Standard_deviations.to_csv(outfilename_Standard_deviations, index = False)
        Times.to_csv(outfilename_Times, index = False)
        Ground_state.to_csv(outfilename_Ground_state, index = False)
        Speedup.to_csv(outfilename_Speedup, index = False)

        files = outfilename_Integrals + " " + outfilename_Errors + " " + outfilename_Standard_deviations \
                + " " + outfilename_Times + " " + outfilename_Ground_state + " " + outfilename_Speedup

        path = "results/monte_carlo/benchmark_test"
        if not os.path.exists(path):
            os.makedirs(path)
        os.system("mv" + " " + files + " " + path)

        print("Integrals:")
        print("------------------------------------------------------------------------------------")
        print(Integrals)
        print("------------------------------------------------------------------------------------")
        print("Ground state")
        print("------------------------------------------------------------------------------------")
        print(Ground_state)
        print("------------------------------------------------------------------------------------")

    if compile == "no":
        outfilename_Integrals = "MC_Integrals.csv"
        outfilename_Errors = "MC_Errors.csv"
        outfilename_Standard_deviations = "MC_Standard_deviations.csv"
        outfilename_Times = "MC_Times.csv"
        outfilename_Ground_state = "MC_Ground_state.csv"
        outfilename_Speedup = "MC_Speedup.csv"


        path ="results/monte_carlo/benchmarks/"

        filenames = [outfilename_Integrals,outfilename_Errors,outfilename_Standard_deviations,outfilename_Times,outfilename_Ground_state,outfilename_Speedup]

        figurename_integrals = "MC_integrals.pdf"
        figurename_STD = "MC_STDs.pdf"
        figurename_times = "MC_timeused.pdf"
        figurename_speedup = "MC_speedup.pdf"
        figurename_Ground_state = "GS_helium.pdf"
        figurename_Errors = "MC_errors.pdf"


        #Plots the results
        exact = 5*np.pi**2/16**2
        integrals = pd.read_csv(path + outfilename_Integrals, header = 0, names = ["N", "BF", "BF_MPI", "IS", "IS_MPI"])
        plt.scatter(np.log10(integrals["N"]), integrals["BF"], label = "Brute force", marker = "x")
        plt.scatter(np.log10(integrals["N"]), integrals["BF_MPI"], label = "Brute force (MPI)", marker = "x")
        plt.scatter(np.log10(integrals["N"]), integrals["IS"], label = "Importance sampling", marker = "x")
        plt.scatter(np.log10(integrals["N"]), integrals["IS_MPI"], label = "Importance sampling (MPI)", marker = "x")
        plt.axhline(y = exact, ls = "--", c = "k", label = "exact")
        plt.legend(fontsize = 12)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        plt.xlabel(r"$\log_{10} N$", fontsize = 16)
        plt.ylabel(r"$I$", fontsize = 16)
        plt.savefig(figurename_integrals, dpi = 1000)
        plt.close()

        Standard_deviations = pd.read_csv(path + outfilename_Standard_deviations, header = 0, names = ["N", "BF", "BF_MPI", "IS", "IS_MPI"])
        plt.plot(np.log10(Standard_deviations["N"]), Standard_deviations["BF"], label = "Brute force")
        plt.plot(np.log10(Standard_deviations["N"]), Standard_deviations["BF_MPI"], label = "Brute force (MPI)")
        plt.plot(np.log10(Standard_deviations["N"]), Standard_deviations["IS"], label = "Importance sampling")
        plt.plot(np.log10(Standard_deviations["N"]), Standard_deviations["IS_MPI"], label = "Importance sampling (MPI)")
        plt.legend(fontsize = 12)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        plt.xlabel(r"$\log_{10} N$", fontsize = 16)
        plt.ylabel(r"$\sigma$", fontsize = 16)
        plt.savefig(figurename_STD, dpi = 1000)
        plt.close()

        Errors = pd.read_csv(path + outfilename_Errors, header = 0, names = ["N", "BF", "BF_MPI", "IS", "IS_MPI"])
        plt.plot(np.log10(Errors["N"]), Errors["BF"], label = "Brute force")
        plt.plot(np.log10(Errors["N"]), Errors["BF_MPI"], label = "Brute force (MPI)")
        plt.plot(np.log10(Errors["N"]), Errors["IS"], label = "Importance sampling")
        plt.plot(np.log10(Errors["N"]), Errors["IS_MPI"], label = "Importance sampling (MPI)")
        plt.legend(fontsize = 12)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        plt.xlabel(r"$\log_{10} N$", fontsize = 16)
        plt.ylabel(r"$\epsilon$", fontsize = 16)
        plt.savefig(figurename_Errors, dpi = 1000)
        plt.close()

        Times = pd.read_csv(path + outfilename_Times, header = 0, names = ["N", "BF", "BF_MPI", "IS", "IS_MPI"])
        plt.plot(np.log10(Times["N"]), Times["BF"], label = "Brute force")
        plt.plot(np.log10(Times["N"]), Times["BF_MPI"], label = "Brute force (MPI)")
        plt.plot(np.log10(Times["N"]), Times["IS"], label = "Importance sampling")
        plt.plot(np.log10(Times["N"]), Times["IS_MPI"], label = "Importance sampling (MPI)")
        plt.legend(fontsize = 12)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        plt.xlabel(r"$\log_{10} N$", fontsize = 16)
        plt.ylabel(r"$\Delta t$ [s]", fontsize = 16)
        plt.savefig(figurename_times, dpi = 1000)
        plt.close()


        Speedup = pd.read_csv(path + outfilename_Speedup, header = 0, names = ["N", "BF", "IS"])
        plt.plot(np.log10(Speedup["N"]), Speedup["BF"], label = "Brute force")
        plt.plot(np.log10(Speedup["N"]), Speedup["IS"], label = "Importance sampling")
        plt.legend(fontsize = 12)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        plt.xlabel(r"$\log_{10} N$", fontsize = 16)
        plt.ylabel("Speedup", fontsize = 16)
        plt.savefig(figurename_speedup, dpi = 1000)
        plt.close()

        Ground_state = pd.read_csv(path + outfilename_Ground_state, header = 0, names = ["N", "H", "error"])
        plt.plot(np.log10(Ground_state["N"]), Ground_state["H"])
        plt.axhline(y = -79, ls = "--", c = "k", label = "Ground state energy")
        plt.axhline(y = -75, ls = "--", c = "r", label = "Exact expectation value")
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        plt.xlabel(r"$\log_{10} N$", fontsize = 16)
        plt.ylabel(r"$\langle H \rangle $", fontsize = 16)
        plt.legend(fontsize = 12)
        plt.savefig(figurename_Ground_state, dpi = 1000)
        plt.close()


        figurenames = figurename_integrals + " " + figurename_STD + " " + figurename_times + " " + figurename_speedup + " " + figurename_Ground_state + " " + figurename_Errors
        path = "results/monte_carlo/plots"
        if not os.path.exists(path):
            os.makedirs(path)

        os.system("mv" + " " + figurenames + " " + path)
