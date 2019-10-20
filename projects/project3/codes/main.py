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

if compilation_instruction == "find_lambda":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    print("executing")

    data = str(input("Produce new data? Type yes or no: "))

    path = "results/benchmarks/"
    main_filename = "Rel_error_legendre.txt"

    h = 0.02
    R = [1 + (i*h) for i in range(1,101)]

    if data == "yes":

        n = 25
        #Radial distance
        integration_method = "1"

        with open(main_filename, "w") as outfile:
            for r in R:
                outfilename = "radial_distance_" + str(r) + "_.txt"
                a = -r
                b = r
                print("executing for r = %f" % r)
                os.system("./main.exe" + " " + outfilename + " "  + integration_method + " " + str(n) +" "+ str(a) +" "+ str(b))
                with open(outfilename, "r") as infile:
                    line = infile.readline()
                    line = line.split()
                    outfile.write(line[2] + "\n")
                os.system("rm"+ " " + outfilename)

        os.system("mv" + " " + main_filename + " " + path)

    rel_error = []

    with open(path +  main_filename, "r") as infile:
        for i in range(len(R)):
            line = infile.readline()
            rel_error.append(float(line))


    index = sum(np.where(rel_error == np.min(rel_error))[0])
    print("Minimum relative error for R = ", R[index])

    figurename = "finding_lamda25.pdf"


    plt.plot(R, rel_error, label="N = 25")
    plt.axvline(x = 2.7, ls = "--", label="$\lambda$ = 2.7")
    plt.xlabel("$\lambda$", fontsize = 16)
    plt.ylabel("$\\xi(\lambda)$", fontsize = 16)
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    plt.legend(fontsize = 16)
    plt.savefig(figurename)
    plt.close()



    os.system("mv" + " " + figurename + " " + path)
    print("Finished!")

if compilation_instruction == "benchmark_legendre":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    print("executing")


    path = "results/benchmarks/"
    main_filename = "finding_N_legendre.txt"


    N = [11,13,15,17,19,21,23,25]
    integration_method = "1"

    a = -3
    b = -a

    with open(main_filename, "w") as outfile:
        for n in N:
            outfilename = "integral_n_" + str(n) + "_.txt"
            print("executing for n = %f" % n)
            os.system("./main.exe" + " " + outfilename + " "  + integration_method + " " + str(n) +" "+ str(a) +" "+ str(b))
            with open(outfilename, "r") as infile:
                line = infile.readline()
                line = line.split()
                outfile.write(line[1] + " " + line[0] + " " + line[2] +" " + line[3] +" "+  "\n")
            os.system("rm"+ " " + outfilename)


    integral = []
    rel_error = []
    timeused = []

    with open(main_filename, "r") as infile:
        for i in range(len(N)):
            line = infile.readline()
            line = line.split()
            integral.append(float(line[1]))
            rel_error.append(float(line[2]))
            timeused.append(float(line[3]))


    data = {\
            "n" : N,\
            "Integration value" : integral,\
            "Relative error" : rel_error, \
            "Time used" : timeused, \
            }



    dataset = pd.DataFrame(data)
    dataset.to_latex(main_filename, index = False)

    os.system("mv" + " " + main_filename + " " + path)




if compilation_instruction == "benchmark_laguerre":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")

    compile = str(input("Produce new data? Type yes or no: "))

    main_filename_3_dim = "benchmark_dim_3_1.csv"
    main_filename_6_dim = "benchmark_dim_6_1.csv"

    if compile == "yes":


        dimensions = [3,6]
        number_of_integration_points = [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]

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

        #Moves the files to appropriate folder
        path = "results/laguerre";
        if not os.path.exists(path):
            os.makedirs(path)
        os.system("mv" + " " +  main_filename_3_dim + " " + main_filename_6_dim + " " + path)

    path = "results/laguerre/"

    #Reading the files generated above

    data_3_dim = pd.read_csv(path + main_filename_3_dim, header = 0, names = ["n_3", "int_val_3", "rel_err_3", "time_3"])
    data_6_dim = pd.read_csv(path + main_filename_6_dim, header = 0, names = ["n_6", "int_val_6", "rel_err_6", "time_6"])

    n = data_3_dim["n_3"]
    exact = 5*np.pi**2/(16*16)

    path = "results/laguerre/"

    data_3_dim = pd.read_csv(path + main_filename_3_dim, header = 0, names = ["n_3", "int_val_3", "rel_err_3", "time_3"])
    data_6_dim = pd.read_csv(path + main_filename_6_dim, header = 0, names = ["n_6", "int_val_6", "rel_err_6", "time_6"])

    n = data_3_dim["n_3"]
    exact = 5*np.pi**2/(16*16)

    #Comparing the results for 3 and 6 dimensions graphically
    figname1 = "integration_value1.pdf"
    figname2 = "time_data1.pdf"
    figname3 = "relative_error1.pdf"


    plt.scatter(n, data_3_dim["int_val_3"], label="3 dimensions")
    plt.scatter(n, data_6_dim["int_val_6"], label="6 dimensions")
    plt.axhline(y = exact, ls = "--", label="Analytical value")
    plt.xlabel("n", fontsize = 14)
    plt.ylabel("Integration value", fontsize = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 14)
    plt.savefig(figname1)
    plt.close()



    plt.plot(np.log10(n), np.log(data_3_dim["time_3"]), label="3 dimensions")
    plt.plot(np.log10(n), np.log(data_6_dim["time_6"]), label="6 dimensions")
    plt.xlabel("$log_{10}(n)$", fontsize = 14)
    plt.ylabel("$log_{10}(Time)$", fontsize = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 14)
    plt.savefig(figname2)
    plt.close()


    plt.plot(n, (data_3_dim["rel_err_3"]), label="3 dimensions")
    plt.plot(n, (data_6_dim["rel_err_6"]), label="6 dimensions")
    plt.xlabel("n", fontsize = 18)
    plt.ylabel("$\epsilon_{relative}$", fontsize = 18)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 14)
    plt.savefig(figname3)
    plt.close()

    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " +  figname1 + " " + figname2 + " " + figname3 + " " + path)


if compilation_instruction == "compare_gauss":
    #Make sure you have already run "benchmark_legendre" and "benckmark_laguerre" before running this section

    #Reading the files from "benchmark_legendre"
    path_legendre = "results/benchmarks/"

    N = [11,13,15,17,19,21,23,25]
    integral_legendre = []
    rel_error_legendre = []
    timeused_legendre = []



    with open(path_legendre + "finding_N_legendre.txt") as infile:
        for j in range(4):
            infile.readline()
        lines = infile.readlines()
        del lines[-1]
        del lines[-1]
        for line in lines:
            line = line.split()
            integral_legendre.append(float(line[2]))
            rel_error_legendre.append(float(line[4]))
            timeused_legendre.append(float(line[6]))

    #Reding the files from "benchmark_laguerre"
    path_laguerre = "results/laguerre/"
    data = pd.read_csv(path_laguerre + "benchmark_dim_3_1.csv" , header = 0, names = ["n", "int_val", "rel_err", "time"])

    integral_laguerre = []
    rel_error_laguerre = []
    time_laguerre = []


    for i in range(1,16,2):
        integral_laguerre.append(data["int_val"][i])
        rel_error_laguerre.append(data["rel_err"][i])
        time_laguerre.append(data["time"][i])

    exact = 5*np.pi**2/16**2

    print(integral_legendre)

    plt.plot(N, integral_legendre, label="Legendre")
    plt.plot(N, integral_laguerre, label="Laguerre")
    plt.axhline(y = exact, ls="--")
    plt.legend()
    plt.show()


    plt.plot(N, timeused_legendre, label="Legendre")
    plt.plot(N, time_laguerre, label="Laguerre")
    plt.legend()
    plt.show()







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

        number_of_monte_carlo_samples = [10**i for i in range(2,7)]
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


        Integrals.to_latex(outfilename_Integrals, encoding='utf-8', escape = False, index = False)
        Errors.to_latex(outfilename_Errors, encoding='utf-8', escape = False, index = False)
        Standard_deviations.to_latex(outfilename_Standard_deviations, encoding='utf-8', escape = False, index = False)
        Times.to_latex(outfilename_Times, encoding='utf-8', escape = False, index = False)
        Ground_state.to_latex(outfilename_Ground_state, encoding='utf-8', escape = False, index = False)
        Speedup.to_latex(outfilename_Speedup, encoding='utf-8', escape = False, index = False)

        files = outfilename_Integrals + " " + outfilename_Errors + " " + outfilename_Standard_deviations \
                + " " + outfilename_Times + " " + outfilename_Ground_state + " " + outfilename_Speedup

        path = "results/monte_carlo/benchmarks"
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
        outfilename_Integrals = "MC_Integrals.txt"
        outfilename_Errors = "MC_Errors.txt"
        outfilename_Standard_deviations = "MC_Standard_deviations.txt"
        outfilename_Times = "MC_Times.txt"
        outfilename_Ground_state = "MC_Ground_state.txt"
        outfilename_Speedup = "MC_Speedup.txt"


        path ="results/monte_carlo/benchmarks/"
        filenames = [outfilename_Integrals,outfilename_Errors,outfilename_Standard_deviations,outfilename_Times,outfilename_Ground_state,outfilename_Speedup]
        for filename in filenames:
            resultater = pd.read_csv(path+filename)
            print(resultater)
