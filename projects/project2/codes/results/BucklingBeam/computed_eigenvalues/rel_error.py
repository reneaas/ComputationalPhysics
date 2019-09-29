import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from linjetilpasning import Linjetilpasning


plt.rc("text", usetex = True)
n = [10, 20, 50, 100, 150, 175, 200, 225, 275, 350]

eig = np.zeros([len(n), 3])


for i in range(len(n)):
    filename = "computed_eigenvalues_BucklingBeam_n_" + str(n[i]) + ".txt"
    with open(filename, "r") as infile:
        lines = infile.readlines()
        eigenvalues = lines[1:4]
        for j in range(len(eigenvalues)):
            eig[i,j] = float(eigenvalues[j].strip("\n"))


real = np.zeros([len(n), 3])

for i in range(len(n)):
    for j in range(2,5):
        real[i,j-2] = (np.pi*j)**2


relative_error = np.zeros([len(n), 3])

for i in range(len(n)):
    for j in range(3):
        relative_error[i,j] = (np.abs(real[i,j] - eig[i,j]))/real[i,j]



rel_err_lambda2 = np.log10(relative_error[:,0])
rel_err_lambda3 = np.log10(relative_error[:,1])
rel_err_lambda4 = np.log10(relative_error[:,2])

nn = [np.log10(N) for N in n]
N = np.linspace(0.5, 3, 1001)

line2 = Linjetilpasning(nn, rel_err_lambda2)
line3 = Linjetilpasning(nn, rel_err_lambda3)
line4 = Linjetilpasning(nn, rel_err_lambda4)


slope_2, I0_2, dslope_2, dI0_2 = line2.linjetilpasning()
slope_3, I0_3, dslope_3, dI0_3 = line3.linjetilpasning()
slope_4, I0_4, dslope_4, dI0_4 = line4.linjetilpasning()

tilpasning2 = slope_2*N + I0_2
tilpasning3 = slope_3*N + I0_3
tilpasning4 = slope_4*N + I0_4

y_error_2 = np.sqrt(dslope_2**2 + dI0_2**2)
y_error_3 = np.sqrt(dslope_3**2 + dI0_3**2)
y_error_4 = np.sqrt(dslope_4**2 + dI0_4**2)


print(slope_2,dslope_2, I0_2, dI0_2)
print(slope_3,dslope_3, I0_3, dI0_3)
print(slope_4,dslope_4, I0_4, dI0_4)



plt.plot(N, tilpasning2,c = "mediumblue" ,label=" $\log_{10}(\epsilon_{rel})$ for $\\lambda_2$")
#plt.errorbar(np.log10(n), rel_err_lambda2, yerr = y_error_2, capsize=6, fmt=".r", label = "Datapoints $\\lambda_2$")
plt.scatter(np.log10(n), rel_err_lambda2, label= "Datapoints $\\lambda_2$")

plt.plot(N, tilpasning3 ,c = "red" ,label=" $\log_{10}(\epsilon_{rel})$ for $\\lambda_3$")
#plt.errorbar(np.log10(n), rel_err_lambda3, yerr = y_error_3, capsize=6, fmt=".y", label = "Datapoints $\\lambda_3$")
plt.scatter(np.log10(n), rel_err_lambda3, label= "Datapoints $\\lambda_3$")

plt.plot(N, tilpasning4, c = "black" ,label=" $\log_{10}(\epsilon_{rel})$ for $\\lambda_2$")
#plt.errorbar(np.log10(n), rel_err_lambda4, yerr = y_error_2, capsize=6, fmt=".k", label = "Datapoints $\\lambda_2$")
plt.scatter(np.log10(n), rel_err_lambda4,c = "grey" ,label= "Datapoints $\\lambda_4$")

plt.xlabel("$log_{10}(n)$" , fontsize = 18)
plt.ylabel("$log_{10}(\epsilon_{rel}$)", fontsize = 18)
plt.legend(fontsize = 14)
plt.xticks(size = 16)
plt.yticks(size = 16)
plt.show()
#plt.savefig("rel_error_bb.png", dpi = 1000)
