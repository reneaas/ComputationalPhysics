import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import os
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
plt.rc("text", usetex = True)

part = str(input("Which part of the project to run: [b,c,d,e,flags] \n"))

if part == "b":
    path = "results/2x2/"
    infilename_exp_ordered = "Expectation_values_n_2_ordered.txt"
    infilename_error_ordered = "Relative_error_n_2_ordered.txt"
    infilename_exp_random =  "Expectation_values_n_2_random.txt"
    infilename_error_random = "Relative_error_n_2_random.txt"

    E_o = []
    E_squared_o = []
    Mabs_o = []
    Mabs_squared_o = []
    M_o = []
    M_squared_o = []
    C_v_o = []
    Chi_o = []


    E_r = []
    E_squared_r = []
    Mabs_r = []
    Mabs_squared_r = []
    M_r = []
    M_squared_r = []
    C_v_r = []
    Chi_r = []


    E_error_o = []
    E_squared_error_o = []
    Mabs_error_o = []
    Mabs_squared_error_o = []
    M_squared_error_o = []
    C_v_error_o = []
    Chi_error_o = []

    E_error_r = []
    E_squared_error_r = []
    Mabs_error_r = []
    Mabs_squared_error_r = []
    M_squared_error_r = []
    C_v_error_r = []
    Chi_error_r = []

    time = []

    with open(path + infilename_exp_ordered, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            time.append(float(values[0]))
            E_o.append(float(values[1]))
            E_squared_o.append(float(values[2]))
            Mabs_o.append(float(values[3]))
            Mabs_squared_o.append(float(values[4]))
            M_o.append(float(values[5]))
            M_squared_o.append(float(values[6]))
            C_v_o.append(float(values[7]))
            Chi_o.append(float(values[8]))

    with open(path + infilename_exp_random, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            E_r.append(float(values[1]))
            E_squared_r.append(float(values[2]))
            Mabs_r.append(float(values[3]))
            Mabs_squared_r.append(float(values[4]))
            M_r.append(float(values[5]))
            M_squared_r.append(float(values[6]))
            C_v_r.append(float(values[7]))
            Chi_r.append(float(values[8]))

    with open(path + infilename_error_ordered, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            E_error_o.append(float(values[1]))
            E_squared_error_o.append(float(values[2]))
            Mabs_error_o.append(float(values[3]))
            Mabs_squared_error_o.append(float(values[4]))
            M_squared_error_o.append(float(values[5]))
            C_v_error_o.append(float(values[6]))
            Chi_error_o.append(float(values[7]))

    with open(path + infilename_error_random, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            E_error_r.append(float(values[1]))
            E_squared_error_r.append(float(values[2]))
            Mabs_error_r.append(float(values[3]))
            Mabs_squared_error_r.append(float(values[4]))
            M_squared_error_r.append(float(values[5]))
            C_v_error_r.append(float(values[6]))
            Chi_error_r.append(float(values[7]))

    beta = 1
    Z_a  = 4 * (3 + np.cosh(8*beta))
    E_a = -32*np.sinh(8*beta)/Z_a
    E_squared_a = 256*np.cosh(8*beta)/Z_a;
    Mabs_a =  8*(np.exp(8*beta) + 2)/Z_a;
    Mabs_squared_a = (32*np.exp(8*beta) + 32)/Z_a;
    chi = (Mabs_squared_a - Mabs_a**2)/1
    Cv = (E_squared_a - E_a**2)/1
    Cv = (256*np.cosh(8) - 32**2*(np.sinh(8))**2/Z_a)/Z_a

    E_ordered = np.array(E_o)
    E_random = np.array(E_r)
    Mabs_ordered = np.array(Mabs_o)
    Mabs_random = np.array(Mabs_r)
    chi_ordered = np.array(Chi_o)
    chi_random = np.array(Chi_r)
    Cv_ordered = np.array(C_v_o)
    Cv_random = np.array(C_v_r)

    E_ordered /= 4
    E_random /= 4
    Mabs_ordered /= 4
    Mabs_random /= 4
    chi_ordered /= 4
    chi_random /= 4
    Cv_ordered /= 4
    Cv_random /= 4

    plt.plot(time[:], E_ordered[:], label = " ordered")
    plt.plot(time[:], E_random[:], label = " random")
    plt.axhline(y = E_a/4, ls = ":", color = "k", label = "analytical")
    plt.xlabel(r"$t$ [cycles/$L^2$]", size = 16)
    plt.ylabel(r"$\langle E\rangle / L^2 $", size = 16)
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    plt.legend(fontsize = 16)
    plt.show()

    plt.plot(time[:], Mabs_ordered[:], label = " ordered")
    plt.plot(time[:], Mabs_random[:], label = " random")
    plt.axhline(y = Mabs_a/4, ls = ":", color = "k", label = "analytical")
    plt.xlabel(r"$t$ [cycles/$L^2$]", size = 16)
    plt.ylabel(r"$\langle |M|\rangle / L^2$", size = 16)
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    plt.legend(fontsize = 16)
    plt.show()


if part == "c":
    T = float(input("Give temperature 1 or 2.4: "))
    path = "results/partC/"
    infilename_ordered = "MC_" + str(int(2e8)) + "_n_20_T_" + str(T) + "_ordered_.txt"
    infilename_random = "MC_" + str(int(2e8)) + "_n_20_T_" + str(T) + "_random_.txt"
    E_ordered = []
    M_ordered = []
    acceptance_ordered = []
    E_random = []
    M_random = []
    acceptance_random = []
    time = []

    with open(path + infilename_ordered, "r") as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            time.append(float(values[0]))
            E_ordered.append(float(values[1]))
            M_ordered.append(float(values[2]))
            acceptance_ordered.append(float(values[3]))

    with open(path + infilename_random, "r") as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            E_random.append(float(values[1]))
            M_random.append(float(values[2]))
            acceptance_random.append(float(values[3]))

    time = [i/1000 for i in time]
    plt.plot(time[:], E_ordered[:], label = "ground state initiation")
    plt.plot(time[:], E_random[:], label = "random initiation")
    plt.xlabel("$t$ [$10^3 \\times$ cycles/$L^2$]", size = 14)
    plt.ylabel(r"$\langle E \rangle / L^2$", size = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 12)
    plt.figure()

    plt.plot(time[:500000], M_ordered[:500000], label = "ground state initiation")
    plt.plot(time[:500000], M_random[:500000], label = "random initiation")
    plt.xlabel("$t$ [$10^3 \\times$ cycles/$L^2$]", size = 14)
    plt.ylabel(r"$\langle |M| \rangle / L^2$", size = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 12)
    plt.figure()

    plt.plot(time[:40000], acceptance_random[:40000], label = "Accepted states (random initiation)")
    plt.plot(time[:40000], acceptance_ordered[:40000], label = "Accepted states (ground state initiation)")
    plt.xlabel("$t$ [$10^3 \\times$ cycles/$L^2$]", size = 14)
    plt.ylabel("Accepted spins", size = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 12)
    plt.show()


if part == "d":
    #N = int(12e7)
    energies = []
    L = 20
    temperature = float(input("Give temperature 1 or 2.4"))
    initial_spin_state = str(input("ordered or random: [o/r] \n"))
    if initial_spin_state == "o":
        initial_spin_state = "ordered"
    if initial_spin_state == "r":
        initial_spin_state = "random"

    infilename =  "boltzmann_distribution_MC_" + str(int(4e7)) + "_T_" +  str(temperature) + "_" + initial_spin_state + "_.txt"
    path = "results/partC/"
    with open(path + infilename, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            values = line.split()
            energies.append(float(values[0]))

    MC = len(energies)
    MC_cycles = np.linspace(0,MC,MC)

    plt.hist(energies, density = True)
    plt.xlabel("$E/J$",size = 16)
    plt.ylabel("$P(E)$",size = 16)
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    plt.show()

    plt.plot(MC_cycles[:200000]/L**2, np.array(energies[:200000]))
    plt.ylabel("$E/J$",size = 16)
    plt.xlabel("$t$ [cycles/$L^2$]",size = 16)
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    plt.show()


if part == "flags":
    path = "results/partE/compilerflag/"
    MC_samples = [100, 1000, 10000, 100000, 1000000, 10000000, 100000000]
    flags = ["No_Flag", "No_Flagp_1", "O2", "O3", "Ofast"]

    for flag in flags:
        time = []
        for MC in MC_samples:
            infilename = "MC_" + str(MC) + "_Flag_" + flag + ".txt"
            file_path = path + infilename
            with open(file_path, "r") as infile:
                lines = infile.readlines()
                for line in lines:
                    values = line.split()
                    time.append(float(values[0]))

        if flag == "O2" or flag == "O3" or flag == "Ofast":
            flag = "-"+ flag
        if flag == "No_Flag":
            flag = "No flag"
        if flag == "No_Flagp_1":
            flag = "No flag (unparallelized)"

        plt.scatter(np.log10(MC_samples), np.log10(time), label = flag, marker = "x")
        plt.xlabel("$\log_{10}(N)$", size = 16)
        plt.xticks(size = 16)
        plt.yticks(size = 16)
        plt.ylabel("$\log_{10}(t)$", size = 16)
        plt.legend(fontsize = 14)

    plt.show()

if part == "e":
    #L = int(input("Lattice size L = "))
    p = 8
    time = 100;
    total_time = 1000000;
    path = "results/partE/total_time_" + str(total_time) + "burn_in_time_" + str(time) + "/"
    my_ranks = [i for i in range(p)]
    Lattice_sizes = [40, 60, 80, 100]
    T_C = []                #To store maximum values.

    fig1 = plt.figure(); figurename_energy = "energies.pdf"
    fig2 = plt.figure(); figurename_magnetization = "magnetization.pdf"
    fig3 = plt.figure(); figurename_chi = "chi.pdf"
    fig4 = plt.figure(); figurename_heat_capacity = "heat_capacity.pdf"
    fig5 = plt.figure(); figurename_interpolation = "interpolated_heat_capacity.pdf"
    fig6 = plt.figure(); figurename_critical_temp = "critical_temp.pdf"
    fig7 = plt.figure(); figurename_std_E = "STD_E.pdf"

    ax1 = fig1.add_subplot(111);
    ax2 = fig2.add_subplot(111);
    ax3 = fig3.add_subplot(111);
    ax4 = fig4.add_subplot(111);
    ax5 = fig5.add_subplot(111);
    ax6 = fig6.add_subplot(111);
    ax7 = fig7.add_subplot(111);

    for L in Lattice_sizes:
        T = []
        E = []
        M = []
        chi = []
        Cv = []
        VarE = []
        infilename = "observables_L_" + str(L) + ".txt";
        file_path = path + infilename
        with open(file_path, "r") as infile:
            lines = infile.readlines()
            for line in lines:
                values = line.split()
                T.append(float(values[0]))
                E.append(float(values[1]))
                M.append(float(values[2]))
                chi.append(float(values[3]))
                Cv.append(float(values[4]))
                VarE.append(float(values[4])*float(values[0])**2*L*L)

        ax1.scatter(T, E, label = str(L) + " x " + str(L), marker = ".")
        ax1.set_xlabel(r"$T$" , size = 14)
        ax1.set_ylabel(r"$\langle E \rangle/L^2J$", size = 14)
        ax1.tick_params(labelsize = 15)
        ax1.legend(fontsize = 12)



        ax2.scatter(T, M, label = str(L) + " x " + str(L), marker = ".")
        ax2.set_xlabel(r"$T$" , size = 14)
        ax2.set_ylabel(r"$\langle |M| \rangle/L^2$", size = 14)
        ax2.tick_params(labelsize = 15)
        ax2.legend(fontsize = 12)


        ax3.scatter(T, chi, label = str(L) + " x " + str(L), marker = ".")
        ax3.set_xlabel(r"$T$" , size = 14)
        ax3.set_ylabel(r"$\chi/L^2$", size = 14)
        ax3.tick_params(labelsize = 15)
        ax3.legend(fontsize = 12)


        ax4.scatter(T, Cv, label = str(L) + " x " + str(L), marker = ".")
        ax4.set_xlabel(r"$T$" , size = 14)
        ax4.set_ylabel(r"$ C_V/L^2$", size = 14)
        ax4.tick_params(labelsize = 15)
        ax4.legend(fontsize = 12)


        #Interpolate dataset and create smoother plots:
        Cs = UnivariateSpline(T,Cv, s=3)
        Ts = np.linspace(2.0,2.399,100)
        Cv_spline = Cs(Ts)
        ax5.plot(Ts, Cv_spline, label = str(L) + " x " + str(L))
        ax5.set_xlabel(r"$T$" , size = 14)
        ax5.set_ylabel(r"$ C_V/L^2$", size = 14)
        ax5.tick_params(labelsize = 15)
        ax5.legend(fontsize = 12)

        index = np.where(np.array(Cv_spline) == max(Cv_spline))
        T_C.append(Ts[np.sum(index)])


        #Plot the standard deviation of E.
        STD_E = [np.sqrt(i) for i in VarE]
        ax7.plot(T,STD_E, label = str(L) + " x " + str(L))
        ax7.set_xlabel(r"$T$" , size = 14)
        ax7.set_ylabel(r"$\sigma_E$", size = 14)
        ax7.tick_params(labelsize = 15)
        ax7.legend(fontsize = 12)


    inverse_L = [1/float(i) for i in Lattice_sizes]
    ax6.scatter(inverse_L, T_C, marker = "x", color = "r", label = "Datapoints")
    ax6.set_xlabel(r"$1/L$", size = 16)
    ax6.set_ylabel(r"$T_C(L)$", size = 16)
    ax6.tick_params(labelsize = 15)
    #ax6.legend(fontsize = 12)

    linear_func = lambda x,a,b: a + b*x             #Form of the function to fit.
    popt, pcov = curve_fit(linear_func, inverse_L, T_C)
    X = np.array(inverse_L)
    L = np.linspace(35,130, 1001)
    X = [1/i for i in L]
    X = np.array(X)
    Y = popt[1]*X + popt[0]
    ax6.plot(X,Y,"-k")
    equation = r"$T_C(L) = 1.152/L + \underbrace{2.266}_{=T_C(\infty)}$"
    ax6.text(0.0080,2.270, equation, {"color":"k", "fontsize": 18})
    ax6.legend(fontsize = 12)
    plt.show()



    fig1.savefig(figurename_energy)
    fig2.savefig(figurename_magnetization)
    fig3.savefig(figurename_chi)
    fig4.savefig(figurename_heat_capacity)
    fig5.savefig(figurename_interpolation)
    fig6.savefig(figurename_critical_temp)
    fig7.savefig(figurename_std_E)

    figurenames = figurename_energy + " " + figurename_magnetization + " " + figurename_chi + " " + figurename_heat_capacity + " " + figurename_interpolation\
                    + " " + figurename_critical_temp + " " + figurename_std_E

    os.system("mv" + " " + figurenames + " " + path)
