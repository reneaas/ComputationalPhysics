import numpy as np
import matplotlib.pyplot as plt
from os import system
f = lambda x: 1 - (1-np.exp(-10))*x - np.exp(-10*x)

Eps = lambda h: M1*h/3 + ((4e-15)/h**2)*M2

zeta = 0.23
M1 = 10**3

M2 = f(zeta)


h = np.linspace(0.000000001, 0.1, 1000001)

eps = Eps(h)

index = np.where(np.log(eps) == min(np.log(eps)))
h_log = np.log10(h)
print(h_log[index])

plt.plot(np.log10(h), np.log(eps), "-r", label="Upper-bound error")
plt.xlabel("log10(h)")
plt.ylabel("log10(Upper-bound-error)")
plt.axvline(x = -5.591563573, label="optimal stepsize h")
plt.legend()

figurename = "upper_bound_error.png"
directory = "~/Documents/skole/comphys/projects/project1/codes/results"
plt.savefig(figurename)
system("mv" + " " + figurename + " " + directory)
