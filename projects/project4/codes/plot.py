import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import os
plt.rc("text", usetex = True)

part = str(input("Which part of the project to run: [b,c,d,e] \n"))

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

    print("E", E_a/4)
    print("M", Mabs_a/4)
    print("chi", chi/4)
    print("Cv", Cv/4)

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


    plt.plot(time[:], E_o[:], label = " ordered")
    plt.plot(time[:], E_r[:], label = " random")
    plt.axhline(y = E_a/4, ls = ":", color = "k", label = "analytical")
    plt.xlabel(r"$t$ [cycles/$L^2$]", size = 16)
    plt.ylabel(r"$\langle E\rangle / L^2 $", size = 16)
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    plt.legend(fontsize = 16)
    #plt.savefig('results/2x2/E.pdf')
    plt.show()
    """
    plt.plot(time[:], M_o[:], label = " ordered")
    plt.plot(time[:], M_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\langle M\rangle$/spins")
    plt.legend()
    plt.savefig('results/2x2/M.pdf')
    plt.close()

    plt.plot(time[:], E_squared_o[:], label = " ordered")
    plt.plot(time[:], E_squared_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\langle E^2\rangle$/spins")
    plt.legend()
    plt.savefig('results/2x2/E_sq.pdf')
    plt.close()

    plt.plot(time[:], M_squared_o[:], label = " ordered")
    plt.plot(time[:], M_squared_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\langle M^2\rangle$/spins")
    plt.legend()
    plt.savefig('results/2x2/M_sq.pdf')
    plt.close()
    """


    plt.plot(time[:], Mabs_o[:], label = " ordered")
    plt.plot(time[:], Mabs_r[:], label = " random")
    plt.axhline(y = Mabs_a/4, ls = ":", color = "k", label = "analytical")
    plt.xlabel(r"$t$ [cycles/$L^2$]", size = 16)
    plt.ylabel(r"$\langle |M|\rangle / L^2$", size = 16)
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    plt.legend(fontsize = 16)
    #plt.savefig('results/2x2/Mabs.pdf')
    plt.show()

    """
    plt.plot(time[:], Mabs_squared_o[:], label = " ordered")
    plt.plot(time[:], Mabs_squared_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\langle |M|^2\rangle$/spins")
    plt.legend()
    plt.savefig('results/2x2/Mabs_sq.pdf')
    plt.close()

    plt.plot(time[:],C_v_o[:], label = " ordered")
    plt.plot(time[:],C_v_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$C_v$/spins")
    plt.legend()
    plt.savefig('results/2x2/Cv.pdf')
    plt.close()

    plt.plot(time[:],Chi_o[:], label = " ordered")
    plt.plot(time[:],Chi_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\chi$/spins")
    plt.legend()
    plt.savefig('results/2x2/Chi.pdf')
    plt.close()


    plt.plot(time[:],E_error_o[:], label = " ordered")
    plt.plot(time[:],E_error_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\epsilon_{\langle E\rangle}$/spins")
    plt.legend()
    plt.savefig('results/2x2/E_error.pdf')
    plt.close()

    plt.plot(time[:],E_squared_error_o[:], label = " ordered")
    plt.plot(time[:],E_squared_error_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\epsilon_{\langle E^2\rangle}$/spins")
    plt.legend()
    plt.savefig('results/2x2/E_sq_error.pdf')
    plt.close()

    plt.plot(time[:],M_squared_error_o[:], label = " ordered")
    plt.plot(time[:],M_squared_error_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\epsilon_{\langle M^2\rangle}$/spins")
    plt.legend()
    plt.savefig('results/2x2/M_sq_error.pdf')
    plt.close()

    plt.plot(time[:],Mabs_error_o[:], label = " ordered")
    plt.plot(time[:],Mabs_error_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\epsilon_{\langle |M|\rangle}$/spins")
    plt.legend()
    plt.savefig('results/2x2/Mabs_error.pdf')
    plt.close()

    plt.plot(time[:],Mabs_squared_error_o[:], label = " ordered")
    plt.plot(time[:],Mabs_squared_error_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\epsilon_{\langle |M|^2\rangle}$/spins")
    plt.legend()
    plt.savefig('results/2x2/Mabs_sq_error.pdf')
    plt.close()

    plt.plot(time[:],C_v_error_o[:], label = " ordered")
    plt.plot(time[:],C_v_error_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\epsilon_{C_v}$/spins")
    plt.legend()
    plt.savefig('results/2x2/Cv_error.pdf')
    plt.close()

    plt.plot(time[:],Chi_error_o[:], label = " ordered")
    plt.plot(time[:],Chi_error_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\epsilon_{\chi}$/spins")
    plt.legend()
    plt.savefig('results/2x2/Chi_error.pdf')
    plt.close()
    """

    values = [100-1, 1000-1, 10000-1, 100000-1, 999999]
    table_values_o = {"t":[], "E":[], "M":[], "heat_cap":[], "chi":[]}
    table_values_r = {"t":[], "E":[], "M":[], "heat_cap":[], "chi":[]}
    for i in values:
        table_values_o["t"].append(np.log10(time[i]))
        table_values_o["E"].append(E_ordered[i])
        table_values_o["M"].append(Mabs_ordered[i])
        table_values_o["heat_cap"].append(Cv_ordered[i])
        table_values_o["chi"].append(chi_ordered[i])
        table_values_r["t"].append(np.log10(time[i]))
        table_values_r["E"].append(E_random[i])
        table_values_r["M"].append(Mabs_random[i])
        table_values_r["heat_cap"].append(Cv_random[i])
        table_values_r["chi"].append(chi_random[i])


    dataset1 = pd.DataFrame(table_values_o)
    dataset1.to_latex("table_2x2_o.tex", index = False, escape = False, encoding = "utf-8")

    dataset2 = pd.DataFrame(table_values_r)
    dataset2.to_latex("table_2x2_r.tex", index = False, escape = False, encoding = "utf-8")



if part == "c":
    T = float(input("Give temperature: "))
    path = "results/partC/"
    infilename_ordered = "MC_" + str(int(4e8)) + "_n_20_T_" + str(T) + "_ordered_.txt"
    infilename_random = "MC_" + str(int(4e8)) + "_n_20_T_" + str(T) + "_random_.txt"
    E_ordered = []
    M_ordered = []
    acceptance_ordered = []
    E_random = []
    M_random = []
    acceptance_random = []
    time = []

    with open(path + infilename_ordered, "r") as infile:
        infile.readline()
        lines = infile.readlinn_spins = n*n;es()
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
    plt.plot(time[:500000], E_ordered[:500000], label = "ground state initiation")
    plt.plot(time[:500000], E_random[:500000], label = "random initiation")
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
    temperature = float(input("Temperature = "))
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

    print(len(energies))
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

if part == "e_old":
    #L = int(input("Lattice size L = "))
    p = 8
    path = "results/partE/"
    my_ranks = [i for i in range(p)]
    Lattice_sizes = [40, 60, 80, 100]

    fig1 = plt.figure(); figurename_energy = "energies.pdf"
    fig2 = plt.figure(); figurename_magnetization = "magnetization.pdf"
    fig3 = plt.figure(); figurename_chi = "chi.pdf"
    fig4 = plt.figure(); figurename_heat_capacity = "heat_capacity.pdf"
    ax1 = fig1.add_subplot(111);
    ax2 = fig2.add_subplot(111);
    ax3 = fig3.add_subplot(111);
    ax4 = fig4.add_subplot(111);

    for L in Lattice_sizes:
        T = []
        E = []
        M = []
        chi = []
        Cv = []
        for my_rank in my_ranks:
            infilename = "observables_my_rank_" + str(my_rank) + "_L_" + str(L) + ".txt"
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

        ax1.plot(T, E, label = str(L) + " x " + str(L))
        ax1.set_xlabel(r"$k_BT$" , size = 14)
        ax1.set_ylabel(r"$\langle E \rangle/J$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax1.legend(fontsize = 12)



        ax2.plot(T, M, label = str(L) + " x " + str(L))
        ax2.set_xlabel(r"$k_BT$" , size = 14)
        ax2.set_ylabel(r"$\langle M \rangle$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax2.legend(fontsize = 12)


        ax3.plot(T, chi, label = str(L) + " x " + str(L))
        ax3.set_xlabel(r"$k_BT$" , size = 14)
        ax3.set_ylabel(r"$\chi$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax3.legend(fontsize = 12)


        ax4.plot(T, Cv, label = str(L) + " x " + str(L))
        ax4.set_xlabel(r"$k_BT$" , size = 14)
        ax4.set_ylabel(r"$ C_V$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax4.legend(fontsize = 12)

    fig1.savefig(figurename_energy)
    fig2.savefig(figurename_magnetization)
    fig3.savefig(figurename_chi)
    fig4.savefig(figurename_heat_capacity)

    figurenames = figurename_energy + " " + figurename_magnetization + " " + figurename_chi + " " + figurename_heat_capacity
    figurepath = "results/partE/plots"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + figurenames + " " + figurepath)

if part == "e_old":
    #L = int(input("Lattice size L = "))
    p = 8
    path = "results/partE/"
    my_ranks = [i for i in range(p)]
    Lattice_sizes = [40, 60, 80, 100]

    fig1 = plt.figure(); figurename_energy = "energies.pdf"
    fig2 = plt.figure(); figurename_magnetization = "magnetization.pdf"
    fig3 = plt.figure(); figurename_chi = "chi.pdf"
    fig4 = plt.figure(); figurename_heat_capacity = "heat_capacity.pdf"
    ax1 = fig1.add_subplot(111);
    ax2 = fig2.add_subplot(111);
    ax3 = fig3.add_subplot(111);
    ax4 = fig4.add_subplot(111);

    for L in Lattice_sizes:
        T = []
        E = []
        M = []
        chi = []
        Cv = []
        for my_rank in my_ranks:
            infilename = "observables_my_rank_" + str(my_rank) + "_L_" + str(L) + ".txt"
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

        ax1.plot(T, E, label = str(L) + " x " + str(L))
        ax1.set_xlabel(r"$k_BT$" , size = 14)
        ax1.set_ylabel(r"$\langle E \rangle/J$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax1.legend(fontsize = 12)



        ax2.plot(T, M, label = str(L) + " x " + str(L))
        ax2.set_xlabel(r"$k_BT$" , size = 14)
        ax2.set_ylabel(r"$\langle M \rangle$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax2.legend(fontsize = 12)


        ax3.plot(T, chi, label = str(L) + " x " + str(L))
        ax3.set_xlabel(r"$k_BT$" , size = 14)
        ax3.set_ylabel(r"$\chi$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax3.legend(fontsize = 12)


        ax4.plot(T, Cv, label = str(L) + " x " + str(L))
        ax4.set_xlabel(r"$k_BT$" , size = 14)
        ax4.set_ylabel(r"$ C_V$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax4.legend(fontsize = 12)

    fig1.savefig(figurename_energy)
    fig2.savefig(figurename_magnetization)
    fig3.savefig(figurename_chi)
    fig4.savefig(figurename_heat_capacity)

    figurenames = figurename_energy + " " + figurename_magnetization + " " + figurename_chi + " " + figurename_heat_capacity
    figurepath = "results/partE/plots"
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("mv" + " " + figurenames + " " + figurepath)

if part == "e":
    #L = int(input("Lattice size L = "))
    p = 8
    time = 100;
    total_time = 1000000;
    path = "results/partE/total_time_" + str(total_time) + "burn_in_time_" + str(time) + "/"
    my_ranks = [i for i in range(p)]
    Lattice_sizes = [40, 60, 80, 100]

    fig1 = plt.figure(); figurename_energy = "energies.pdf"
    fig2 = plt.figure(); figurename_magnetization = "magnetization.pdf"
    fig3 = plt.figure(); figurename_chi = "chi.pdf"
    fig4 = plt.figure(); figurename_heat_capacity = "heat_capacity.pdf"
    ax1 = fig1.add_subplot(111);
    ax2 = fig2.add_subplot(111);
    ax3 = fig3.add_subplot(111);
    ax4 = fig4.add_subplot(111);

    for L in Lattice_sizes:
        T = []
        E = []
        M = []
        chi = []
        Cv = []
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

        ax1.scatter(T, E, label = str(L) + " x " + str(L), marker = ".")
        ax1.set_xlabel(r"$k_BT$" , size = 14)
        ax1.set_ylabel(r"$\langle E \rangle/L^2J$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax1.legend(fontsize = 12)



        ax2.scatter(T, M, label = str(L) + " x " + str(L), marker = ".")
        ax2.set_xlabel(r"$k_BT$" , size = 14)
        ax2.set_ylabel(r"$\langle |M| \rangle/L^2$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax2.legend(fontsize = 12)


        ax3.scatter(T, chi, label = str(L) + " x " + str(L), marker = ".")
        ax3.set_xlabel(r"$k_BT$" , size = 14)
        ax3.set_ylabel(r"$\chi/L^2$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax3.legend(fontsize = 12)


        ax4.scatter(T, Cv, label = str(L) + " x " + str(L), marker = ".")
        ax4.set_xlabel(r"$k_BT$" , size = 14)
        ax4.set_ylabel(r"$ C_V/L^2$", size = 14)
        plt.xticks(size = 14)
        plt.yticks(size = 14)
        ax4.legend(fontsize = 12)

    fig1.savefig(figurename_energy)
    fig2.savefig(figurename_magnetization)
    fig3.savefig(figurename_chi)
    fig4.savefig(figurename_heat_capacity)

    figurenames = figurename_energy + " " + figurename_magnetization + " " + figurename_chi + " " + figurename_heat_capacity

    os.system("mv" + " " + figurenames + " " + path)
