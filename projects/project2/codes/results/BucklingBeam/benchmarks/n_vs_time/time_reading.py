import matplotlib.pyplot as plt
import numpy as np
import os
from linjetilpasning import Linjetilpasning
import seaborn as sns

plt.rc("text", usetex = True)
n = [10, 20, 50, 100, 150]

Jacobi_time = []
Arma_time = []

for i in range(5):
    filename = "bb_time_" + str(n[i]) + ".txt"
    infile = open(filename)
    for line in infile:
        time = line.split()
        Jacobi_time.append(float(time[1]))
        Arma_time.append(float(time[3]))
    #os.system("rm" + " " + filename)

n = [np.log10(N) for N in n]

Jacobi_time = [np.log10(i) for i in Jacobi_time]
Arma_time = [np.log10(i) for i in Arma_time]

Jacobi = np.array(Jacobi_time)
Arma = np.array(Arma_time)

Line_jacobi = Linjetilpasning(n, Jacobi)
Line_arma = Linjetilpasning(n, Arma)

slope_jacobi, I0_jacobi, dslope_jacobi, dI0_jacobi = Line_jacobi.linjetilpasning()
slope_arma, I0_arma, dslope_arma, dI0_arma = Line_arma.linjetilpasning()


N = np.linspace(0.5, 2.5, 1001)
time_jacobi = slope_jacobi*N + I0_jacobi
time_arma = slope_arma*N + I0_arma

y_error_jacobi = np.sqrt(dslope_jacobi**2 + dI0_jacobi**2)
y_error_arma = np.sqrt(dslope_arma**2 + dI0_arma**2)



plt.plot(N, time_jacobi,c = "mediumblue" ,label=" $\log_{10}(t)$ (Jacobi) ")
plt.plot(N, time_arma,c = "k" ,label="$\log_{10}(t) $ (Armadillo)")
plt.errorbar(n, Jacobi, yerr = y_error_jacobi, capsize=4, fmt=".r", label = "Datapoints Jacobi")
plt.errorbar(n, Arma, yerr = y_error_arma, capsize=4,fmt=".b", label = "Datapoints Armadillo")
plt.xlabel(r"$\log_{10} (n)$" ,fontsize = 18)
plt.ylabel(r"$\log_{10} (t)$", fontsize = 18)
plt.legend(fontsize = 14)
plt.xticks(size = 16)
plt.yticks(size = 16)
plt.show()
#plt.savefig("time_vs_n.pdf", dpi = 1000)
