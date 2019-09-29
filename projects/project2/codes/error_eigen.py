import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import seaborn
plt.rc("text", usetex=True)


rho_maxes = [2.0,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,\
                5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,8.0,9.0,10.0,11.0,\
                12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,25.0,30.0,35.0,40.0,45.0,50.0,100.0,200.0]

eig_1 = []
eig_2 = []
eig_3 = []
eig_4 = []




for rho_max in rho_maxes:  #Picks out the first four eigenvalues for each value of rho and adds them to lists

    path = "results/QM_OneElectron/computed_eigenvalues/"
    filename = "computed_eigenvalues_QM_OneElectron_350_" + str(float(rho_max)) + ".txt"
    outfilename = path + filename
    with open(outfilename,"r") as outfile:
        lines = outfile.readlines()
        line1 = float(lines[1])
        eig_1.append(line1)
        line2 = float(lines[2])
        eig_2.append(line2)
        line3 = float(lines[3])
        eig_3.append(line3)
        line4 = float(lines[4])
        eig_4.append(line4)

#Turning the lists into arrays
eig_1 = np.array(eig_1)
eig_2 = np.array(eig_2)
eig_3 = np.array(eig_3)
eig_4 = np.array(eig_4)

#The analytical eigenvalues
eig_ana = np.array([3, 7, 11, 15])


#Calculates the relative error and finds the smallest error
error_1 = np.abs((eig_1-eig_ana[0])/eig_ana[0])
min_1 = np.where(error_1==np.min(error_1))
print(min_1)
error_2 = np.abs((eig_2-eig_ana[1])/eig_ana[1])
min_2 = np.where(error_2==np.min(error_2))
print(min_2)
error_3 = np.abs((eig_3-eig_ana[2])/eig_ana[2])
min_3 = np.where(error_3==np.min(error_3))
print(min_3)
error_4 = np.abs((eig_4-eig_ana[3])/eig_ana[3])
min_4 = np.where(error_4==np.min(error_4))
print(min_4)


tot_err = error_1 + error_2 + error_3 + error_4
min_tot = np.where(tot_err==np.min(tot_err))
print(min_tot)



plt.plot(rho_maxes, error_1, label = "$\lambda_1$")
plt.plot(rho_maxes, error_2, label = "$\lambda_2$")
plt.plot(rho_maxes, error_3, label = "$\lambda_3$")
plt.plot(rho_maxes, error_4, label = "$\lambda_4$")
plt.plot(rho_maxes[11], error_1[11],"o", label = "Optimal $\\rho_{max}$ for $\lambda_1$ = %.2f"%rho_maxes[11])
plt.plot(rho_maxes[15], error_2[15],"o",label = "Optimal $\\rho_{max}$ for $\lambda_2$ = %.2f"%rho_maxes[15])
plt.plot(rho_maxes[19], error_3[19],"o",label = "Optimal $\\rho_{max}$ for $\lambda_3$ = %.2f"%rho_maxes[19])
plt.plot(rho_maxes[23], error_4[23],"o",label = "Optimal $\\rho_{max}$ for $\lambda_4$ = %.2f"%rho_maxes[23])
plt.legend(fontsize = 14)
plt.xticks(size = 16)
plt.yticks(size = 16)
plt.xlabel("$\\rho_{max}$", size = 16)
plt.ylabel("$\epsilon_{rel}$ for $\lambda_i$, $i = 1,2,3,4$", size = 16)
plt.show()


plt.plot(rho_maxes[9:25], error_1[9:25]/(10**-5),"-", label = "$\lambda_1$")
plt.plot(rho_maxes[9:25], error_2[9:25]/(10**-5),"-",label = "$\lambda_2$")
plt.plot(rho_maxes[9:25], error_3[9:25]/(10**-5),"-",label = "$\lambda_3$")
plt.plot(rho_maxes[9:25], error_4[9:25]/(10**-5),"-",label = "$\lambda_4$")
plt.plot(rho_maxes[11], error_1[11]/(10**-5),"o", label = "Optimal $\\rho_{max}$ for $\lambda_1$ = %.2f"%rho_maxes[11])
plt.plot(rho_maxes[15], error_2[15]/(10**-5),"o",label = "Optimal $\\rho_{max}$ for $\lambda_2$ = %.2f"%rho_maxes[15])
plt.plot(rho_maxes[19], error_3[19]/(10**-5),"o",label = "Optimal $\\rho_{max}$ for $\lambda_3$ = %.2f"%rho_maxes[19])
plt.plot(rho_maxes[23], error_4[23]/(10**-5),"o",label = "Optimal $\\rho_{max}$ for $\lambda_4$ = %.2f"%rho_maxes[23])
plt.legend(fontsize = 14)
plt.xticks(size = 16)
plt.yticks(size = 16)
plt.xlabel("$\\rho_{max}$", size = 16)
plt.ylabel("$\epsilon_{rel}/10^{-4}$ for $\lambda_i$, $i = 1,2,3,4$", size = 16)
plt.show()



plt.plot(rho_maxes,tot_err)
plt.plot(rho_maxes[23], tot_err[23],"o",label = "Optimal $\\rho_{max}$ for least total error = %.2f"%rho_maxes[23])
plt.legend(fontsize = 14)
plt.xticks(size = 16)
plt.yticks(size = 16)
plt.xlabel("$\\rho_{max}$", size = 16)
plt.ylabel("Total $\epsilon_{rel}$", size=16)
plt.show()


plt.plot(rho_maxes[21:25], tot_err[21:25]/(10**-4),"-")
plt.plot(rho_maxes[23], tot_err[23]/(10**-4),"o",label = "Optimal $\\rho_{max}$ for least total error = %.2f"%rho_maxes[23])
plt.legend(fontsize = 14)
plt.xticks(size = 16)
plt.yticks(size = 16)
plt.xlabel("$\\rho_{max}$", size = 16)
plt.ylabel("Total $\epsilon_{rel}/10^{-4}$", size=16)
plt.show()
