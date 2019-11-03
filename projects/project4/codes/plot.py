import numpy as np
import matplotlib.pyplot as plt
import sys
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



    plt.plot(time[:], E_o[:], label = " ordered")
    plt.plot(time[:], E_r[:], label = " random")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\langle E\rangle$/spins")
    plt.legend()
    plt.savefig('results/2x2/E.pdf')
    plt.close()

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

    plt.plot(time[:], Mabs_o[:], label = " ordered")
    plt.plot(time[:], Mabs_r[:], label = " ordered")
    plt.xlabel(r"$t$ [cycles/spins]")
    plt.ylabel(r"$\langle |M|\rangle$/spins")
    plt.legend()
    plt.savefig('results/2x2/Mabs.pdf')
    plt.close()

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


if part == "c":
    T = float(input("Give temperature: "))
    path = "results/partC/"
    infilename_ordered = "MC_" + str(int(4e6)) + "_n_20_T_" + str(T) + "_ordered_.txt"
    infilename_random = "MC_" + str(int(4e6)) + "_n_20_T_" + str(T) + "_random_.txt"
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


    plt.plot(time[:], E_ordered[:], label = "ground state initiation")
    plt.plot(time[:], E_random[:], label = "random initiation")
    plt.xlabel(r"$t$ [cycles/spins]", size = 14)
    plt.ylabel(r"$\langle E \rangle $/spins", size = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 12)
    plt.figure()

    plt.plot(time[:], M_ordered[:], label = "ground state initiation")
    plt.plot(time[:], M_random[:], label = "random initiation")
    plt.xlabel("$t$ [cycles/spins]", size = 14)
    plt.ylabel(r"$\langle |M| \rangle $/spins", size = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 12)
    plt.figure()

    plt.plot(time[:], acceptance_random[:], label = "Accepted states (random initiation)")
    plt.plot(time[:], acceptance_ordered[:], label = "Accepted states (ground state initiation)")
    plt.xlabel("$t$ [cycles/spins]", size = 14)
    plt.ylabel("Accepted spins", size = 14)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.legend(fontsize = 12)
    plt.show()

if part == "d":
    N = int(4e6)
    energies = []
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

    MC = len(energies)
    MC_cycles = np.linspace(0,MC,MC)

    #plt.hist(energies, 400 + 1, density = True)
    #plt.show()
    plt.plot(MC_cycles[N:], energies[N:])
    plt.show()

if part == "e":
    L = int(input("Lattice size L = "))
    my_ranks = [i for i in range(8)]
    T = []
    E = []
    M = []
    chi = []
    Cv = []


    for my_rank in my_ranks:
        infilename = "observables_my_rank_" + str(my_rank) + "_L_" + str(L) + ".txt"
        with open(infilename, "r") as infile:
            lines = infile.readlines()
            for line in lines:
                values = line.split()
                T.append(float(values[0]))
                E.append(float(values[1]))
                M.append(float(values[2]))
                chi.append(float(values[3]))
                Cv.append(float(values[4]))


    plt.plot(T, E, label = str(L) + " x " + str(L))
    plt.xlabel(r"$k_BT$" )
    plt.ylabel(r"$\langle E \rangle/J$")
    plt.show()

    plt.plot(T, M, label = str(L) + " x " + str(L))
    plt.xlabel(r"$k_BT$" )
    plt.ylabel(r"$\langle M \rangle$")
    plt.show()

    plt.plot(T, chi, label = str(L) + " x " + str(L))
    plt.xlabel(r"$k_BT$" )
    plt.ylabel(r"$\chi$")
    plt.show()

    plt.plot(T, Cv, label = str(L) + " x " + str(L))
    plt.xlabel(r"$k_BT$" )
    plt.ylabel(r"$ C_V$")
    plt.show()
